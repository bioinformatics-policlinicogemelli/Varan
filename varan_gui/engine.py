"""The container engine (Docker, or Podman): find it, start it, run Varan.

Everything goes through the engine's own command line program, so the app
needs no Python package to talk to it.
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path

# No console window flashing up for every command on Windows.
NO_WINDOW = 0x08000000 if sys.platform == "win32" else 0

CANDIDATES = {
    "docker": {
        "win32": [r"C:\Program Files\Docker\Docker\resources\bin\docker.exe"],
        "darwin": ["/usr/local/bin/docker", "/opt/homebrew/bin/docker",
                   "/Applications/Docker.app/Contents/Resources/bin/docker",
                   str(Path.home() / ".docker" / "bin" / "docker")],
        "linux": ["/usr/bin/docker", "/usr/local/bin/docker"],
    },
    "podman": {
        "win32": [r"C:\Program Files\RedHat\Podman\podman.exe"],
        "darwin": ["/opt/podman/bin/podman", "/opt/homebrew/bin/podman",
                   "/usr/local/bin/podman"],
        "linux": ["/usr/bin/podman", "/usr/local/bin/podman"],
    },
}

DOWNLOAD_PAGE = {
    "docker": "https://www.docker.com/products/docker-desktop/",
    "podman": "https://podman-desktop.io/",
}


def _platform() -> str:
    return "win32" if sys.platform == "win32" else "darwin" if sys.platform == "darwin" else "linux"


def find_cli(kind: str) -> str | None:
    found = shutil.which(kind)
    if found:
        return found
    for candidate in CANDIDATES.get(kind, {}).get(_platform(), []):
        if Path(candidate).is_file():
            return candidate
    return None


class Engine:
    def __init__(self, kind: str = "docker") -> None:
        self.kind = kind if kind in CANDIDATES else "docker"
        self.cli = find_cli(self.kind)

    # An app started from the macOS Finder gets a minimal PATH, without the
    # folder of the engine's helpers (docker-credential-desktop, needed even
    # to pull a public image).
    def _env(self) -> dict[str, str]:
        env = dict(os.environ)
        extra = [str(Path(self.cli).parent)] if self.cli else []
        if _platform() == "darwin":
            extra += ["/usr/local/bin", "/opt/homebrew/bin",
                      "/Applications/Docker.app/Contents/Resources/bin"]
        env["PATH"] = os.pathsep.join([*extra, env.get("PATH", "")])
        return env

    def call(self, *args: str, timeout: float = 60) -> subprocess.CompletedProcess:
        if not self.cli:
            raise FileNotFoundError(f"{self.kind} is not installed")
        return subprocess.run(
            [self.cli, *args], capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=timeout, env=self._env(),
            stdin=subprocess.DEVNULL, creationflags=NO_WINDOW, check=False)

    @property
    def installed(self) -> bool:
        return self.cli is not None

    def running(self) -> bool:
        if not self.cli:
            return False
        try:
            return self.call("info", "--format", "{{.ServerVersion}}", timeout=30).returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            return False

    def start(self) -> bool:
        """Start Docker Desktop / the Podman machine; False when the app
        can't do it (Linux services need an administrator)."""
        system = _platform()
        try:
            if self.kind == "podman" and system != "linux":
                subprocess.Popen([self.cli, "machine", "start"], env=self._env(),
                                 stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL, creationflags=NO_WINDOW)
                return True
            if self.kind == "docker" and system == "win32":
                root = os.environ.get("ProgramFiles", r"C:\Program Files")
                exe = Path(root) / "Docker" / "Docker" / "Docker Desktop.exe"
                if exe.is_file():
                    subprocess.Popen([str(exe)], creationflags=NO_WINDOW)
                    return True
                return False
            if self.kind == "docker" and system == "darwin":
                return subprocess.run(["open", "-a", "Docker"], check=False).returncode == 0
        except OSError:
            return False
        return False

    def wait(self, timeout: float = 240, cancelled: Callable[[], bool] = lambda: False) -> bool:
        end = time.monotonic() + timeout
        while time.monotonic() < end and not cancelled():
            if self.running():
                return True
            time.sleep(3)
        return False

    def image_id(self, image: str) -> str:
        """Local id of `image` ("" if it isn't there)."""
        try:
            done = self.call("image", "inspect", "--format", "{{.Id}}", image)
        except (OSError, subprocess.TimeoutExpired):
            return ""
        return done.stdout.strip() if done.returncode == 0 else ""

    def platform_args(self) -> list[str]:
        """The image is built for x86-64: Apple silicon Macs run it emulated."""
        if _platform() == "darwin" and platform.machine() == "arm64":
            return ["--platform", "linux/amd64"]
        return []

    def stream(self, args: list[str], on_line: Callable[[str], None],
               on_start: Callable[[subprocess.Popen], None] = lambda p: None) -> int:
        """Run an engine command, passing each output line to `on_line`;
        returns its exit code."""
        proc = subprocess.Popen(
            [self.cli, *args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL, env=self._env(), creationflags=NO_WINDOW,
            text=True, encoding="utf-8", errors="replace", bufsize=1)
        on_start(proc)
        assert proc.stdout is not None
        for line in proc.stdout:
            on_line(line.rstrip("\r\n"))
        return proc.wait()

    def pull(self, image: str, on_line: Callable[[str], None]) -> int:
        return self.stream(["pull", *self.platform_args(), image], on_line)

    def kill(self, name: str) -> None:
        try:
            self.call("kill", name, timeout=30)
        except (OSError, subprocess.TimeoutExpired):
            pass
