# Varan: proposte di nuove analisi e funzioni

Raccolta delle proposte emerse da sei profili di utente (bioinformatico, oncologo MTB, biologa wet-lab,
direttore di laboratorio, ricercatore clinico, revisore di rivista), unificate e senza doppioni.
Lista da discutere una per una. Le proposte approvate vengono implementate sul branch `dag2-features`.

Effort: **S** piccolo (giorni), **M** medio (1-3 settimane), **L** grande (mesi).
Stato: `da discutere` / `approvata` / `implementata` / `scartata` / `rimandata`.

Principio comune: Varan non rifà ciò che cBioPortal fa già (OncoPrint, Kaplan-Meier, group comparison,
lollipop); produce i dati e gli attributi che rendono utili quelle funzioni, più ciò che il portale non fa.

---

# DA VALUTARE

## A. Sicurezza dell'interpretazione (tumor-only e biopsia liquida)

### 2. Flag "sospetto germinale" — S — `rimandata` (richiede una discussione a parte, più complessa)
Variante P/LP in ClinVar + gene della lista ESMO (BRCA1, BRCA2, PALB2, MLH1, MSH2, MSH6, RET sempre; altri per
tipo tumorale) + VAF compatibile con germinale (soglia configurabile) → suggerisce consulenza genetica.
Oggi `VAF_EXCLUDE_BANDS` ([SigMA]) toglie queste varianti invece di segnalarle.
Output: colonna MAF, attributo campione, case list "germline follow-up suggested".
Rif.: Kuzbari 2023 Ann Oncol doi:10.1016/j.annonc.2022.12.003.
Rischi: purezza alta o LOH ingannano la VAF; è un sospetto, non una diagnosi; consenso informato.

### 3. "Negativo non vuol dire negativo" — M — `rimandata` (analisi in "Discussione #3" in fondo; servono file MetricsOutput reali)
Per ogni gene chiave (ESCAT I): coperto / non coperto / sotto soglia (MetricsOutput + BED del pannello);
flag su purezza bassa, CNV non valutabili, frazione tumorale ctDNA bassa; frase "assenza non conclusiva".
Rif.: Pascual 2022 Ann Oncol (ESMO ctDNA) doi:10.1016/j.annonc.2022.05.520.

## B. Azionabilità e Molecular Tumor Board

### 8. Report pre-MTB per paziente — M — `da discutere`
Una pagina HTML/PDF: qualità, biomarcatori, alterazioni per livello, VUS separate, flag (CH, germinale),
spazio per la raccomandazione. Etichetta Research Use Only.
Rif.: ESMO PMWG MTB recommendations, Ann Oncol 2025 doi:10.1016/j.annonc.2025.02.009.

### 9. Ricerca di trial clinici aperti — L — `da discutere`
API ClinicalTrials.gov v2 / EU CTIS; escono solo gene, alterazione, tumore. Opzionale (richiede rete).
Rif.: Klein 2022 npj Precis Oncol (MatchMiner) doi:10.1038/s41698-022-00312-5.

### 10. Registro delle decisioni MTB e indicatori — M — `da discutere`
Data discussione, raccomandazione, terapia avviata, motivo del mancato avvio, esito; indicatori ESMO.
Contesto: DM 30/05/2023 sui MTB regionali (da verificare sulla fonte primaria).

## C. Dati clinici e longitudinali

### 11. Dizionario dati clinici validato — S/M — `da discutere`
YAML al posto di HEADER_* in conf.ini; campi standard (OS/PFS, stadio, ECOG, fumo) con tipo e dominio,
controllati in preflight.

### 14. Left truncation e coerenza delle date — S — `da discutere`
`DAYS_DX_TO_SEQ`, avviso nel report, export entry/exit/event, controlli (prelievo dopo diagnosi...).
Rif.: Brown 2022 JAMA Oncol doi:10.1001/jamaoncol.2021.5153.

### 15. Pseudonimizzazione — M — `da discutere`
HMAC con chiave locale, date in giorni relativi, controllo k-anonymity, avviso su ID che sembrano codici
ospedalieri. Prerequisito per condivisione ed export (37).

## D. Analisi genomiche

### 16. Attributi di pathway (Sanchez-Vega) — S/M — `da discutere`
Per campione: `PATHWAY_RTK_RAS` = Altered / Not altered / Not profiled (tiene conto del pannello), solo
alterazioni oncogeniche OncoKB.
Rif.: Sanchez-Vega 2018 Cell doi:10.1016/j.cell.2018.03.035.

### 17. TMB calcolato da Varan e armonizzato — M — `da discutere`
TMB dal MAF / Mb del BED, con IC di Poisson, accanto a quello del vendor; calibrazione opzionale verso WES.
Rif.: Merino 2020 JITC doi:10.1136/jitc-2019-000147; Vega 2021 Ann Oncol doi:10.1016/j.annonc.2021.09.016.

### 18. Controlli di coerenza tra biomarcatori — S/M — `da discutere`
MSI-H senza firma MMR, SBS3 discordante con BRCA/GIS → "da rivalutare". Dipende da SigMA.
Rif.: Gulhan 2019 Nat Genet doi:10.1038/s41588-019-0390-2.

### 19. Purezza, LOH, zigosità — M (SGZ leggero) / L (PureCN) — `da discutere`
Rif.: Sun 2018 PLoS Comput Biol doi:10.1371/journal.pcbi.1005965; Riester 2016 doi:10.1186/s13029-016-0060-z.

### 20. Co-occorrenza/esclusività con DISCOVER — M — `da discutere`
Rif.: Canisius 2016 Genome Biol doi:10.1186/s13059-016-1114-x.

### 21. Controllo identità / scambio campioni (somalier) — M — `da discutere`
Richiede BAM. Rif.: Pedersen 2020 Genome Med doi:10.1186/s13073-020-00761-2.

## E. Qualità e governo del laboratorio

### 22. Metriche QC come attributi + soglie — S/M — `da discutere`
Copertura mediana, % esoni >=100X, contaminazione, PASS/FAIL da MetricsOutput; segnala o esclude.
Rif.: Ewels 2016 (MultiQC) doi:10.1093/bioinformatics/btw354.

### 23. Archivio storico e cruscotto KPI — M — `da discutere`
TAT, tasso fallimenti, carte di controllo, deriva TMB/MSI tra lotti.
Rif.: Tawari 2017 (ChronQC) doi:10.1093/bioinformatics/btx843; Jennings 2017 J Mol Diagn doi:10.1016/j.jmoldx.2017.01.011.

### 24. Fascicolo ISO 15189 + validazione su materiali di riferimento — S/M — `da discutere`
Manifest con checksum e versioni; `--validate-against <truth set>` → sensibilità/PPV.
Rif.: Roy 2018 J Mol Diagn doi:10.1016/j.jmoldx.2017.11.003.

## F. Prodotti per ricerca e pubblicazione

### 25. Table 1 + diagramma di flusso STROBE — S/M — `da discutere`
Celle con n<5 oscurate. Rif.: von Elm 2007 STROBE PMID 17947786.

### 26. Methods automatico + figure di coorte — S/M — `da discutere`
Dalla provenance; oncoplot, distribuzioni TMB/MSI, heatmap pathway in SVG/PDF.
Rif.: Mayakonda 2018 (maftools) doi:10.1101/gr.239244.118.

### 27. Confronto con GENIE/TCGA (geni in comune) — M/L — `da discutere`
Rif.: AACR GENIE 2017 Cancer Discov doi:10.1158/2159-8290.CD-17-0151.

### 28. Query di fattibilità — S/M — `da discutere`
`varan feasibility "KRAS G12C AND LUAD"` → n, IC95%, potenza indicativa.

## G. Facilità d'uso

### 29. Interfaccia grafica — M (Gooey) / L (Streamlit) — `da discutere`
Quattro pulsanti, preflight immediato, errori in linguaggio semplice.
Rif.: Ustjanzew 2021 (cbpManager) doi:10.1186/s12911-021-01719-z.

### 30. Profili di configurazione + template Excel con OncoTree — S — `da discutere`
"TSO500 solido", "TSO500 liquid", "Guardant"; menu a tendina; suggerimenti ("intendevi LUAD?").

## H. Infrastruttura e pubblicazione (Varan 2)

### 31. Test pytest + CI GitHub Actions — S/M — `da discutere`

### 32. Immagine pubblicata, bioconda, versioni fissate (anche ClinVar) — M — `da discutere`

### 33. GRCh38 — L — `da discutere`
Rif.: Hanssen 2024 (nf-core/sarek 3) doi:10.1093/nargab/lqae031.

### 34. Benchmark SEQC2 + concordanza con referti vendor su >=100 campioni — M/L — `da discutere`
Rif.: Fang 2021 Nat Biotechnol doi:10.1038/s41587-021-00993-6 (da ricontrollare); Deveson 2021 doi:10.1038/s41587-021-00857-z.

### 35. Altri vendor (Oncomine, Archer, FoundationOne, DRAGEN, nf-core/sarek) — M/L — `da discutere`

### 36. Genome Nexus come annotatore opzionale — M — `da discutere`
Rif.: de Bruijn 2022 JCO CCI doi:10.1200/CCI.21.00144.

## I. Interoperabilità

### 37. Export AACR GENIE — M — `da discutere` (dipende da 15)

### 38. RO-Crate — S — `da discutere`
Rif.: Soiland-Reyes 2022 Data Science doi:10.3233/DS-210053.

### 39. FHIR Genomics / mCODE, GA4GH VRS, Phenopackets, Beacon — L — `da discutere`
Rif.: Osterman 2020 JCO CCI doi:10.1200/CCI.20.00059; Wagner 2021 Cell Genomics doi:10.1016/j.xgen.2021.100027.

### 40. Connettore di ingresso dal LIMS (CSV) — S/M — `da discutere`

## In sospeso da proposte già implementate
- #13 decisione 5: "Clearance" (oggi "Not detected" separato da Decrease). Da decidere insieme alla #3.
- #6: manca la tabella ESCAT vera del laboratorio, da convertire nel formato di escat_template.tsv.

## Nota regolatoria trasversale
Le funzioni 3, 6, 8, 9, se usate per scegliere la terapia del singolo paziente, rendono plausibile la
qualifica di software IVD (IVDR, probabile classe C). Etichetta Research Use Only su ogni output,
versione di ogni base di conoscenza nel report, validazione umana, verifica licenza OncoKB per uso clinico.

## Ordine suggerito
(Aggiornato togliendo le proposte già fatte: 1, 4, 5, 6, 7, 12, 13.)
1. Rapide (S): 2, 16, 22, 26, 30.
2. Studio clinico: 11, 14, 15, 25.
3. Solidità per Varan 2: 31, 32, 17, 34.
4. Dopo: 29, 8, 33, 35.

---

## Discussione #3 — "Negativo non vuol dire negativo"

### Il problema
Se in un campione non compare una mutazione in KRAS, oggi in cBioPortal quel campione risulta "non alterato" in
KRAS. Ma "non trovata" può voler dire cose diverse:
1. il gene non è nel pannello;
2. il gene è nel pannello ma quella zona è stata letta male (copertura bassa);
3. il campione ha troppo poco tumore (purezza o frazione di ctDNA bassa): una mutazione vera avrebbe una VAF
   sotto il limite di rilevazione del test;
4. quel tipo di dato è fallito il QC (per esempio TMB o MSI non valutabili, CNV non valutabili nel plasma);
5. il tipo di dato non è stato proprio testato (per esempio niente RNA, quindi niente fusioni).

### Cosa Varan copre già
- Punto 1: gene panel matrix (`[GenePanel]`, anche per tipo con `[GenePanel.Plasma]`) → cBioPortal mostra
  "not profiled" per i geni fuori pannello.
- Punto 5: colonne `cnv_tested` / `fusion_tested` e case list dei campioni profilati.
- Non c'è niente sui punti 2, 3 e 4. Varan oggi **non legge MetricsOutput.tsv** (verificato nel codice).

### Cosa si potrebbe aggiungere: quattro livelli
| Livello | Cosa | Input | Effort | Affidabilità |
|---|---|---|---|---|
| A | QC di campione per tipo di dato (PASS/FAIL dal vendor) | MetricsOutput.tsv TSO500; equivalente Guardant da capire | S-M | alta: è il QC ufficiale del vendor |
| B | Sensibilità attesa da purezza / frazione tumorale | colonna TC già in sample.tsv; LOD del test in conf.ini | S | media: stima grezza (clone eterozigote, nessuna correzione per copy number) |
| C | Copertura dei geni/hotspot chiave | gVCF o file di copertura per esone, se esistono nei vostri output | M-L | alta se c'è il dato |
| D | Frase nel report per campione ("assenza non conclusiva per: ...") | output di A-C | S | quanto i livelli sotto |

Il livello A coincide in parte con la proposta #22 (metriche QC come attributi): farli insieme.

### Domande aperte
1. Quali livelli: A + B + D subito, C dopo?
2. QC FAIL di un tipo di dato: solo segnalare, o anche togliere il campione dai "profilati" per quel tipo?
3. Come arriva MetricsOutput: colonna opzionale in sample.tsv o ricerca automatica nella cartella del run?
4. Soglie QC: quelle di Illumina come default in conf.ini, sovrascrivibili per tipo?
5. LOD per tipo di campione: valori di default?
6. Geni chiave: tutti quelli del pannello, una lista in conf.ini, o una lista curata per tumore?
7. Guardant: cosa danno come QC e come frazione tumorale? Servono file di esempio.
8. Avete gVCF o file di copertura per esone dai run TSO500? Servono file di esempio per il livello C.

---

# FATTE

## A. Sicurezza dell'interpretazione (tumor-only e biopsia liquida)

### 1. Segnalazione di possibile ematopoiesi clonale (CH) — S/M — `implementata` (dag2-features, commit 720f4d0)
Segnalare (mai rimuovere di default) le varianti che potrebbero venire dai globuli bianchi e non dal tumore.
Chi la chiede: 5 profili su 6. Vedi "Discussione #1" in fondo.
Rif.: Razavi 2019 Nat Med doi:10.1038/s41591-019-0652-7; Hu 2018 Clin Cancer Res doi:10.1158/1078-0432.CCR-18-0143;
Bolton 2020 Nat Genet doi:10.1038/s41588-020-00710-0.

## B. Azionabilità e Molecular Tumor Board

### 4. Attributi di azionabilità per campione + case list derivate — S — `implementata` (dag2-features, commit 49ac20d; vedi "Discussione #4-#7" in fondo)
`HIGHEST_ONCOKB_LEVEL`, `N_ACTIONABLE`; case list TMB-high, MSI-H, OncoKB livello 1, tessuto/plasma.
Riepilogo di coorte nel report ("il 34% ha un'alterazione di livello 1-2"). Salvare versione/data OncoKB.
Rif.: Chakravarty 2017 JCO PO doi:10.1200/PO.17.00011.

### 5. Riquadro biomarcatori tumore-agnostici — S — `implementata` (dag2-features, commit 49ac20d; vedi "Discussione #4-#7" in fondo)
NTRK1-3, RET, ALK/ROS1, BRAF V600E, FGFR, ERBB2, MSI-H, TMB-H, HRD/GIS, con stato "non valutabile".
Rif.: Mosele 2024 Ann Oncol doi:10.1016/j.annonc.2024.04.005; Marabelle 2020 Lancet Oncol doi:10.1016/S1470-2045(20)30445-9.

### 6. Livelli ESCAT e AMP/ASCO/CAP accanto a OncoKB — M — `implementata` (dag2-features, commit 49ac20d; vedi "Discussione #4-#7" in fondo)
Tabella curata gene/alterazione x tumore → ESCAT, datata e citata. Mappatura OncoTree → tumore ESMO.
Rif.: Mateo 2018 Ann Oncol doi:10.1093/annonc/mdy263; Li 2017 J Mol Diagn doi:10.1016/j.jmoldx.2016.10.002.
Rischi: curatela continua, responsabilità sul contenuto.

### 7. CIViC come base di conoscenza aperta — M — `implementata` (dag2-features, commit 49ac20d; vedi "Discussione #4-#7" in fondo)
Dump CC0 offline, alternativa senza token a OncoKB; evidenziare discordanze OncoKB/CIViC.
Rif.: Griffith 2017 Nat Genet doi:10.1038/ng.3774; Horak 2022 Genet Med doi:10.1016/j.gim.2022.01.001.

## C. Dati clinici e longitudinali

### 12. Timeline (trattamenti, prelievi, stato/risposta) — M — `implementata` (dag2-features, commit 870c998; OS/PFS esclusi; vedi "Discussione #12" in fondo)
`data_timeline_*` in giorni relativi a un'origine; il validatore incluso supporta già il formato.
update/extract/remove devono propagarli.
Rif.: Lavery 2022 JCO CCI (GENIE BPC) doi:10.1200/CCI.21.00105; RECIST 1.1 doi:10.1016/j.ejca.2008.10.026.

### 13. ctDNA seriale — M/L — `implementata` (dag2-features, commit 3c26558; decisione 6 in 870c998; decisione 5 in sospeso)
VAF max/media per prelievo, variazioni, varianti emergenti di resistenza, clearance.
Rif.: Vega 2022 JCO PO (ctMoniTR) doi:10.1200/PO.21.00372.

---

## Discussione #1 — Possibile ematopoiesi clonale (CH)

### Cosa dice la letteratura
- **Razavi et al. 2019, Nat Med** (MSK; doi:10.1038/s41591-019-0652-7). 124 pazienti metastatici con biopsia
  tissutale contemporanea + 47 controlli sani; cfDNA e globuli bianchi (WBC) appaiati sequenziati su 508 geni,
  ~2 Mb, >60.000x. Risultati: il 53,2% delle mutazioni nel cfDNA dei pazienti (81,6% nei controlli) aveva
  caratteristiche di CH; la maggior parte delle mutazioni con VAF <1% era presente nei WBC; geni più
  frequenti DNMT3A, TET2, PPM1D, TP53. Conclusione degli autori: serve il sequenziamento appaiato cfDNA-WBC
  per interpretare correttamente le varianti.
- **Non è un metodo gold standard da implementare**: è lo studio di riferimento che dimostra il problema.
  Il gold standard è il **sequenziamento dei WBC dello stesso prelievo**: una variante presente nei WBC è CH.
  Senza WBC si può solo stimare.
- **Hu et al. 2018, Clin Cancer Res** (DFCI; doi:10.1158/1078-0432.CCR-18-0143): falsi positivi nel plasma da CH
  in NSCLC. La maggior parte delle JAK2 V617F, parte delle TP53 e rare KRAS G12X nel cfDNA venivano dai WBC.
- **Bolton et al. 2020, Nat Genet** (doi:10.1038/s41588-020-00710-0): CH in ~24.000 pazienti oncologici; le
  CH nei geni del danno al DNA (PPM1D, TP53, CHEK2) aumentano dopo chemio/radioterapia. Fonte per la lista
  geni/hotspot.
- **Metodi senza WBC**: MetaCH (npj Precis Oncol 2025, doi:10.1038/s41698-025-00921-w; classificatore ML
  addestrato anche sui dati Razavi); Fairchild 2023 Sci Transl Med doi:10.1126/scitranslmed.abm8729
  (classificazione da solo cfDNA su 4324 campioni). Sono modelli di ricerca, non standard clinici.
- **Guardant360** applica già un proprio filtro CHIP proprietario (dalla scheda tecnica): sugli input Guardant
  una parte delle CH potrebbe essere già tolta a monte. Da verificare cosa resta nei file che riceviamo.

### Perché "VAF bassa" da sola non funziona
- Nel liquido tutte le VAF sono basse; la VAF di una variante tumorale dipende dalla frazione di ctDNA,
  quella di una CH no (dipende dalla dimensione del clone nel sangue).
- Una soglia fissa in conf.ini è arbitraria: con ctDNA basso anche i driver tumorali sono allo 0,5%.
- Normalizzare sul campione (VAF della variante / VAF massima delle varianti tumorali sicure) aiuta ma non
  decide: un rapporto basso può essere CH **o** un subclone tumorale; un rapporto alto non esclude CH se il
  ctDNA è quasi assente (la CH può essere la variante a VAF più alta del campione).
- Conclusione: la VAF è un modificatore, non il criterio principale.

### Proposta di implementazione (livelli di evidenza, dal più forte)
1. **WBC appaiato disponibile** (nuova colonna opzionale in sample.tsv con VCF dei WBC): variante presente nei
   WBC → `CH confermata`. Unico caso in cui si può "chiamare" con sicurezza.
2. **Evidenza interna allo studio** (Varan ha già più campioni per paziente):
   - stessa variante presente nel **tessuto** dello stesso paziente → `tumorale (presente nel tessuto)`;
   - assente nel tessuto con posizione coperta → indizio a favore di CH;
   - **plasma seriale**: VAF stabile mentre le varianti tumorali cambiano con la terapia → indizio a favore di CH.
3. **Gene e tipo di variante** (lista versionata `ch_genes.tsv` con fonte):
   - geni tipicamente CH (DNMT3A, TET2, ASXL1, PPM1D tronc. esone 6, JAK2 V617F, SRSF2, SF3B1, U2AF1, CBL,
     GNB1...): indizio forte;
   - geni "doppi" (TP53, ATM, CHEK2, KRAS, IDH1/2, GNAS): indizio solo se combinato con 2 o VAF relativa bassa.
4. **VAF relativa al campione**: modificatore, soglie in conf.ini con default documentati.

### Nomi e presentazione nell'output
- Mai "CHIP" o "CH" come diagnosi senza WBC. Colonna MAF `CH_EVIDENCE` con valori:
  `CH confermata (WBC)` / `Sospetta CH` / `Compatibile con CH` / `Tumorale (supportata da tessuto)` / vuoto.
- Colonna `CH_REASON` con il motivo ("gene DNMT3A; assente nel tessuto S12; VAF 0,4% vs max 8%").
- In cBioPortal: custom driver tier nel MAF (`cbp_driver_tiers`, da verificare il nome esatto del campo) →
  filtro nel portale per includere/escludere le sospette CH. Attributo campione `N_CH_SUSPECT` e geni.
- Nel report HTML: sezione per campione con le varianti segnalate e il motivo, più la frase
  "Senza sequenziamento dei globuli bianchi l'origine non può essere confermata".
- Nessuna rimozione di default. Opzione in conf.ini per escludere le `Sospetta CH` dai calcoli a valle
  (TMB, SigMA, pathway, frequenze), dichiarata nel report.
- Default: attivo per campioni Liquid; per il tessuto opzionale (CH anche da infiltrato ematico, VAF bassa).

### Affidabilità
- Con WBC: alta (è il gold standard).
- Con tessuto dello stesso paziente: buona per varianti presenti nel tessuto; "assente nel tessuto" è solo
  indizio (eterogeneità, campionamento).
- Solo gene + VAF: euristica; non sappiamo quantificarla senza validazione. Proposta: validare su un set con
  WBC (se il laboratorio ne ha) o su dati pubblici con annotazione WBC, prima di dichiarare numeri.

### Vincoli attuali di Varan (verificati nel codice)
- Il tipo di campione è per studio, non per campione: `[Sample_Type] TYPE` = Solid o Liquid. Tessuto e
  plasma dello stesso paziente stanno nello stesso studio solo dopo un merge con `-u`. Per usare il tessuto
  come evidenza (livello 2) serve una colonna `SAMPLE_TYPE` per campione (tissue/plasma), conservata nel merge.
- L'adapter Guardant non legge oggi nessun campo legato a CH. Da controllare sui file Guardant reali
  (VCF FILTER/INFO, XML) se Guardant marca le varianti CH: se sì, le usiamo come evidenza.

### Tipo di campione per singolo campione (decisione: colonna opzionale con override + warning)

Dove Varan oggi usa `[Sample_Type] TYPE` (verificato nel codice):
| Punto | File | Cosa cambia tra Solid e Liquid |
|---|---|---|
| Lettura MSI dal CombinedVariantOutput | tsv.py `get_msi_tmb` | sezione diversa del file |
| Classificazione MSI | walk_clinical.py | Liquid: `THRESHOLD_MSI_LIQUID`; Solid: `THRESHOLD_MSI` |
| Vendor | varan.py `reconcile_sample_type` | PIPELINE=guardant implica Liquid |
| Filtro VAF | filter_clinvar.py | `t_VAF_min` unico (il commento in conf.ini dice 0.03 solid / 0.001 liquid) |
| SigMA | sigma_filter.py | parametri unici, split per tipo già indicato come domanda aperta |
| Report e provenance | write_report.py | `sample_type` unico; nel merge mostrato solo se uguale nei due studi |
| Gene panel | gene_panels.py | un pannello per tipo di dato per tutta la corsa |

Proposta:
1. Colonna opzionale in sample.tsv (nome da decidere, vedi sotto), valori Solid/Liquid. Se compilata vince su
   `[Sample_Type] TYPE` per quel campione, con warning che elenca i campioni sovrascritti. Celle vuote =
   valore dello studio.
2. Scritta sempre in data_clinical_sample.txt per ogni campione (anche quando viene dal valore generale),
   così il tipo sopravvive a update/extract/remove. Per gli studi vecchi senza colonna, nel merge `-u` si
   riempie dal `sample_type` della provenance di quello studio.
3. Confini: ogni sezione di conf.ini che dipende dal tipo accetta una sezione di override per tipo, come
   già si fa con `[Vendor.guardant]`: `[Filters.Liquid]`, `[MSI.Liquid]`, `[TMB.Liquid]`, `[SigMA.Liquid]`,
   `[GenePanel.Liquid]` (e `.Solid`). Per ogni campione Varan legge prima la sezione del suo tipo, poi
   quella generale. Un'unica funzione `setting(section, key, sample_id)` sostituisce le letture dirette.
4. Se una corsa contiene entrambi i tipi e manca la sezione di override per uno dei due, warning: "i
   campioni Liquid usano le soglie generali [Filters]". MSI resta com'è (ha già le due soglie).
5. Coerenza vendor: un campione Guardant marcato Solid → errore (o warning) esplicito.
6. Report: tipo "Mixed (12 Solid, 5 Liquid)" e tabella delle soglie effettive per tipo.

Flusso consigliato: una corsa per tipo (ognuna con la sua conf), poi merge con `-u`. La colonna serve per i
casi misti in una sola corsa e soprattutto per conservare il tipo dopo il merge, che è ciò che serve al
livello 2 della CH.

Nome della colonna: cBioPortal usa già `SAMPLE_TYPE` per Primary/Metastasis/Recurrence (icone nella patient
view). Usarlo per Solid/Liquid creerebbe conflitto con chi lo compila già con quel significato. Proposta:
`SAMPLE_SOURCE` (Solid/Liquid) o `SPECIMEN_TYPE` (Tissue/Plasma).

### Decisioni prese
- [x] Nome colonna tipo di campione: **SPECIMEN_TYPE**
- [x] Valori: **Tissue / Plasma**, anche per i confini (vedi specifica)
- [x] Sezioni di override per tipo: **approvate** (`[Filters.Plasma]` ecc.)
- [x] Guardant marcato Tissue: **errore, la corsa si blocca**
- [x] Categorie in inglese, colonne `CH_EVIDENCE` + `CH_REASON`
- [x] VCF WBC come input opzionale: **sì**
- [x] Altri campioni dello stesso paziente come evidenza: **sì**
- [x] Lista geni: bozza nostra (Bolton 2020 + Razavi 2019 + Hu 2018), nessuna revisione esterna (uso di ricerca)
- [x] Soglie VAF relativa in conf.ini con default documentati, dichiarate euristiche
- [x] Attiva **solo sui campioni Plasma**
- [x] Esclusione dai calcoli a valle: opzione, **spenta di default**
- [x] cBioPortal: custom driver tier + attributi `N_CH_SUSPECT`, `CH_SUSPECT_GENES` (nome campo tier da verificare)
- [x] Avvertenza fissa nel report
- [x] Nessuna validazione con WBC possibile: nessuna percentuale di accuratezza dichiarata

**Stato #1: implementata** sul branch `dag2-features` (commit 720f4d0, non ancora pushato).

### Note di implementazione (differenze rispetto alla specifica)
- **Niente custom driver tier.** Verificato sul validatore di cBioPortal: `cbp_driver_tiers` serve a marcare
  le varianti come *driver* (menu "Mutation Color" dell'OncoPrint); usarlo per la CH le farebbe apparire come
  driver. Al loro posto: colonne nel **namespace "CH"** di cBioPortal (`CH.EVIDENCE`, `CH.REASON`, con
  `namespaces: CH` nel meta file), che il portale mostra nelle tabelle delle mutazioni. Il filtro per campione
  resta con `N_CH_SUSPECT` / `CH_SUSPECT_GENES`. Il validatore accetta tutto.
- **Colonne MAF**: `CH.EVIDENCE` / `CH.REASON` invece di `CH_EVIDENCE` / `CH_REASON` (il punto è richiesto
  dal namespace).
- **VCF WBC**: colonna opzionale `wbc_path` in sample.tsv. Con il WBC presente, una variante assente dal WBC
  non viene segnalata (motivo "not in matched WBC").
- **Esclusione da SigMA** (`EXCLUDE_FROM_DOWNSTREAM = True`): usa solo l'evidenza disponibile per il singolo
  campione (WBC, gene, VAF relativa), perché SigMA gira prima che gli altri campioni siano pronti.
  TMB calcolato e pathway non esistono ancora: si agganceranno quando arrivano (#16, #17).
- **Dopo un merge `-u`** le segnalazioni vengono ricalcolate con tutti i campioni dello studio unito: il tessuto
  arrivato con un batch conferma o smentisce il plasma dell'altro. Le "CH confirmed (WBC)" restano anche se il
  VCF WBC non viene ripassato.
- **Default `SERIAL_FOLD = 2`** per il plasma seriale (VAF della variante che cambia meno di 2 volte mentre la
  VAF tumorale cambia almeno 2 volte).
- conf.ini del repo: `TYPE = Plasma` (prima `Liquid`, stesso significato).

### Test eseguiti
- Regressione completa (10 scenari, anche Snakemake): unica differenza attesa la nuova colonna `SPECIMEN_TYPE`
  (più le colonne CH sullo studio Guardant plasma); validatore cBioPortal OK.
- Test dedicato tessuto + plasma dello stesso paziente, con WBC e plasma seriale: tutte le categorie escono come
  previsto, identiche con varan.py e Snakemake; filtro VAF diverso per plasma via `[Filters.Plasma]`; pannelli
  diversi per tipo via `[GenePanel.Plasma]`.
- Campioni Plasma con `PIPELINE = illumina_solid`: la corsa si blocca con l'errore.
- Merge `-u` di uno studio tissue e uno plasma "vecchio" (senza colonna): `SPECIMEN_TYPE` recuperato dalla
  provenance, CH ricalcolata con il tessuto.

### Specifica approvata

**Tipo di campione (Tissue / Plasma)**
- Vocabolario unico ovunque: `Tissue` e `Plasma`.
  - Colonna opzionale `SPECIMEN_TYPE` in sample.tsv: vince su `[Sample_Type] TYPE` per quel campione,
    con warning che elenca i campioni sovrascritti. Cella vuota = valore dello studio.
  - `[Sample_Type] TYPE` accetta `Tissue` / `Plasma`. I vecchi `Solid` / `Liquid` restano accettati come
    sinonimi (Solid = Tissue, Liquid = Plasma) con warning di deprecazione: i conf.ini esistenti funzionano.
  - Sezioni di override: `[Filters.Tissue]`, `[Filters.Plasma]`, `[TMB.Plasma]`, `[SigMA.Plasma]`,
    `[GenePanel.Plasma]`... Per ogni campione si legge prima la sezione del suo tipo, poi quella generale.
  - `THRESHOLD_MSI_LIQUID` resta valido come alias; il nuovo posto è `[MSI.Plasma] THRESHOLD_MSI`.
  - I `PIPELINE` `illumina_solid` / `illumina_liquid` restano (nomi di prodotto TSO500) e implicano
    Tissue / Plasma; `guardant` implica Plasma.
- Perché Tissue/Plasma: "Solid" è ambiguo (tumore solido vs tessuto: anche il plasma viene da un tumore
  solido); Tissue/Plasma descrive il materiale. Altri materiali futuri (liquor, urine) si aggiungono come
  nuovi valori senza rinominare nulla.
- Scritto sempre in data_clinical_sample.txt per ogni campione. Merge `-u` con studi vecchi senza colonna:
  riempita dal `sample_type` della provenance (Solid → Tissue, Liquid → Plasma).
- Report: "Mixed (12 Tissue, 5 Plasma)" e tabella delle soglie effettive per tipo.
- Corsa mista senza sezione di override per un tipo presente → warning con i valori generali usati.

**Segnalazione CH (solo campioni Plasma)**
- Colonne MAF:
  - `CH_EVIDENCE`: `CH confirmed (WBC)` / `CH suspected` / `CH compatible` / `Tumor (tissue-supported)` / vuoto
  - `CH_REASON`: motivo leggibile, es. "DNMT3A; absent in tissue S12 (covered); VAF 0.4% vs sample max 8%"
- Regole, dalla più forte:
  1. Variante presente nel VCF WBC del campione → `CH confirmed (WBC)`.
     Nuova colonna opzionale in sample.tsv per il VCF WBC (nome da fissare in implementazione).
  2. Variante presente in un campione Tissue dello stesso paziente nello studio → `Tumor (tissue-supported)`.
  3. Gene tier 1 (tipicamente CH) → `CH suspected`; le evidenze aggiuntive (assente nel tessuto coperto,
     VAF stabile nel plasma seriale) arricchiscono il motivo.
  4. Gene tier 2 (doppio) → `CH compatible`; sale a `CH suspected` solo con un'evidenza in più
     (assente nel tessuto coperto, VAF relativa sotto soglia, VAF stabile nel seriale mentre i driver cambiano).
  5. VAF relativa (VAF / VAF massima delle varianti non-CH del campione): solo modificatore.
- conf.ini, nuova sezione `[CH]`:
  - `ENABLED = True` (agisce solo sui campioni Plasma)
  - `GENE_LIST = ch_genes.tsv` (file versionato nel repo, con le fonti)
  - `RELATIVE_VAF_MAX = 0.25` (default prudente: sotto un quarto della VAF tumorale massima)
  - `EXCLUDE_FROM_DOWNSTREAM = False` (se True: escluse da TMB calcolato, SigMA, pathway; dichiarato nel report)
- cBioPortal: custom driver tier nel MAF per filtrare le sospette CH; attributi campione `N_CH_SUSPECT`,
  `CH_SUSPECT_GENES`.
- Report HTML: sezione per campione Plasma con varianti segnalate e motivo; avvertenza fissa:
  "Without sequencing of matched white blood cells the origin of these variants cannot be confirmed.
  Research use only."
- Nessuna percentuale di accuratezza dichiarata.

**Bozza lista geni `ch_genes.tsv`** (fonti: Razavi 2019, Hu 2018, Bolton 2020; da rifinire in implementazione)

| Tier | Gene | Varianti | Nota |
|---|---|---|---|
| 1 | DNMT3A | non-sinonime; hotspot R882 | gene CH più frequente |
| 1 | TET2 | troncanti e missense | |
| 1 | ASXL1 | troncanti | |
| 1 | PPM1D | troncanti esoni 5-6 | arricchito dopo chemio/radio |
| 1 | JAK2 | V617F | Hu 2018: quasi sempre da WBC |
| 1 | SRSF2 | P95 | |
| 1 | SF3B1 | K700E e hotspot | |
| 1 | U2AF1 | S34, Q157 | |
| 1 | ZRSR2 | troncanti | |
| 1 | GNB1 | K57 | |
| 1 | CBL | missense dominio RING | |
| 1 | BCOR, BCORL1, STAG2 | troncanti | |
| 2 | TP53 | tutte | anche tumorale; arricchito dopo terapia |
| 2 | CHEK2 | tutte | anche germinale |
| 2 | ATM | tutte | rischio di falsa indicazione a PARP-inibitori |
| 2 | KRAS, NRAS | hotspot | Hu 2018: rare KRAS G12X da WBC |
| 2 | IDH1, IDH2 | R132 / R140, R172 | anche tumorali (glioma, colangiocarcinoma) |
| 2 | GNAS | R201 | anche tumorale |

I geni assenti dal pannello (TSO500 ctDNA, Guardant360) semplicemente non vengono mai segnalati.

---

## Discussione #13 — ctDNA seriale (più prelievi di plasma dello stesso paziente)

### L'idea
Un paziente fa più biopsie liquide nel tempo (prima della terapia, durante, alla progressione). Confrontandole si
vede se il tumore risponde (il ctDNA scende o sparisce) e se compaiono mutazioni nuove (spesso di resistenza).

### Letteratura
- ctMoniTR (Vega et al. 2022, JCO Precis Oncol, doi:10.1200/PO.21.00372): 200 pazienti NSCLC trattati con
  immunoterapia, 5 studi; la riduzione del ctDNA durante la terapia si associa a OS e PFS migliori.
- Non esiste una definizione condivisa di "risposta molecolare". Le più usate: calo >= 50% della VAF (massima o
  media) rispetto al basale; clearance = nessuna variante più rilevabile (più rara, più specifica).
- ESMO 2022 (Pascual et al.): il monitoraggio ctDNA è ancora ricerca, non pratica clinica.

### Cosa c'è già
- Varan gestisce più campioni per paziente; con la #1 sa quali sono Plasma e segnala le varianti da CH.
- cBioPortal ha già, nella vista paziente, il grafico della VAF per variante tra campioni ("Genomic Evolution"),
  ordinato per data se lo studio ha i prelievi nella timeline. Varan oggi non scrive date né timeline.

### Cosa aggiungerebbe Varan
| Pezzo | Cosa | Effort |
|---|---|---|
| 1. Tempo | giorni del prelievo (o data convertita in giorni) per campione → ordine dei prelievi, file timeline "SPECIMEN" per cBioPortal | S |
| 2. Metriche per prelievo | VAF massima, VAF media, n. varianti, escluse le varianti CH e quelle di aspetto germinale | S |
| 3. Confronto col basale | variazione % rispetto al primo prelievo; categoria: Baseline / Response (calo >= 50%) / Stable / Increase / Not detected | S |
| 4. Varianti per variante | nuova rispetto al prelievo precedente, persistente, non più rilevata (solo se il gene è coperto in entrambi) | M |
| 5. Resistenza | lista curata di mutazioni di resistenza per terapia (es. EGFR C797S, ESR1, KRAS/NRAS dopo anti-EGFR) | L |
| 6. Report | tabella e grafico per paziente | S-M |

### Domande aperte
1. Come arriva il tempo: data assoluta o giorni? (privacy: in uscita solo giorni relativi)
2. Basale: il primo prelievo in ordine di tempo o una colonna esplicita?
3. Metrica: VAF massima, media, o entrambe? Escludere CH e germinali: sì?
4. Soglie: Response >= 50% di calo? Increase? Configurabili in conf.ini?
5. "Clearance": la chiamiamo solo con i dati di QC/LOD (#3), altrimenti "Not detected"?
6. Terapia: serve sapere quando inizia? Collegarla alla #12 (timeline trattamenti) o lasciarla fuori per ora?
7. Resistenza (pezzo 5): ora o dopo?
8. Pannelli diversi tra prelievi: confrontare solo i geni coperti da entrambi?

### Decisioni #13
- [x] 1. Tempo: colonna `COLLECTION_DATE` (data, YYYY-MM-DD o DD/MM/YYYY). Le date restano in
  `service/collection_dates.tsv` nella cartella dello studio; in cBioPortal solo i giorni dal primo prelievo.
- [x] 2. Basale: il primo prelievo di plasma; colonna opzionale `BASELINE` (yes) per sceglierne un altro.
- [x] 3. Metrica configurabile: `[Serial] METRIC` = max (default) o mean; `EXCLUDE_CH`, `GERMLINE_VAF_BANDS`.
- [x] 4. Soglie configurabili con default dalla letteratura, usati se il campo è vuoto: Decrease = calo > 50%,
  Increase = aumento > 20%, Intermediate il resto (ctMoniTR, analisi TKI, Clin Cancer Res 2025;31(11):2162;
  50% anche nell'analisi immunoterapia, Vega 2022).
- [ ] 5. **In sospeso:** "Clearance". Per ora "Not detected" separato da Decrease (ctMoniTR conta la clearance come
  Decrease). Da decidere insieme alla #3 (QC/LOD).
- [x] 6. Terapia (commit 870c998, con la #12): se il paziente ha terapie in treatments.tsv e nessun BASELINE,
  ogni prelievo è confrontato con l'ultimo prelievo di plasma fino all'inizio della terapia in corso ("Not evaluable
  (no pre-treatment sample)" se non c'è; fuori terapia "Baseline" se fa da basale, altrimenti "No treatment
  ongoing"). Nuove colonne RESPONSE_BASELINE, TREATMENT_AT_SAMPLE, DAYS_ON_TREATMENT. Resistenza segnalata solo se
  il paziente ha ricevuto una terapia della classe prima del prelievo; le altre restano solo nel report
  (Received = no). Senza terapie: come prima (Received = unknown). BASELINE esplicito vince sempre.
- [x] 7. Lista di resistenza implementata ora: `resistance_mutations.tsv`, 15 voci, ogni riferimento verificato su
  PubMed (DOI e PMID); solo le varianti nominate negli abstract. Escluse con motivazione scritta nel file: NRAS con
  anti-EGFR, altre mutazioni ectodominio EGFR, ESR1 L536/E380Q, altre ALK, NTRK3, FGFR2 diverse da V564F, AR
  W742C/H875Y, KIT secondarie, PIK3CA. Le reversioni BRCA1/2 non si riconoscono in automatico: segnalata ogni
  variante BRCA nuova con nota "verificare che ripristini il frame". KRAS hotspot e reversioni BRCA solo se nuove.
  Nel report: metodo, soglie, fonti, tabella con DOI.
- [x] 8. Pannelli diversi tra prelievi: "New" solo se il gene era coperto prima, altrimenti
  "New (not covered before)"; "Lost" solo se il gene è coperto nel prelievo successivo.

### Cosa produce (implementato)
- `data_timeline_specimen.txt` → cBioPortal ordina i prelievi nella vista paziente ("Genomic Evolution").
- Attributi campione: DAYS_FROM_FIRST_SAMPLE, CTDNA_MAX_VAF_PCT, CTDNA_MEAN_VAF_PCT, CTDNA_N_VARIANTS,
  CTDNA_CHANGE_PCT, MOLECULAR_RESPONSE, CTDNA_N_NEW, CTDNA_N_LOST, RESISTANCE_VARIANTS.
- MAF, namespace SERIAL: SERIAL.STATUS (Baseline / Persistent / New / New (not covered before)), SERIAL.RESISTANCE.
- Ricalcolato in update/extract/remove: due batch uniti con -u danno lo stesso risultato di una corsa unica
  (testato). Validatore cBioPortal OK, Snakemake identico a varan.py.

---

## Discussione #12 — Timeline clinica (terapie, stato di malattia)

### L'idea
cBioPortal sa mostrare, per ogni paziente, una linea del tempo con terapie, risposte, progressioni e prelievi
(file "timeline"). Con la #13 Varan scrive già i **prelievi** (SPECIMEN). Mancano gli eventi clinici: quando inizia e
finisce ogni terapia, con quale farmaco e in quale linea, quando c'è progressione o risposta.

### Perché serve
- Leggere il ctDNA seriale nel contesto: un calo dopo l'inizio della terapia è una risposta; prima no.
- Decisione 6 della #13: con la terapia nota, il basale diventa "l'ultimo prelievo prima dell'inizio terapia" e una
  mutazione di resistenza viene segnalata solo se il paziente ha ricevuto quella terapia.
- È la base di ogni analisi "alterazione → risposta al farmaco" (stesso schema di AACR GENIE BPC; Lavery 2022,
  JCO CCI, doi:10.1200/CCI.21.00105) e della futura #11 (OS/PFS).

### Cosa aggiungerebbe Varan
| Pezzo | Cosa | Effort |
|---|---|---|
| 1. Input | file opzionale `treatments.tsv` (una riga per evento) passato con -i dopo patient.tsv | S |
| 2. Date → giorni | come per i prelievi: date in ingresso, salvate in `service/`, in cBioPortal solo giorni | S |
| 3. Timeline TREATMENT | `data_timeline_treatment.txt` (AGENT, linea, classe di terapia) | S |
| 4. Timeline STATUS | `data_timeline_status.txt` (progressione, risposta RECIST, decesso) | S |
| 5. Collegamento con #13 | basale rispetto all'inizio terapia; resistenza filtrata sulla terapia ricevuta | M |
| 6. Controlli di coerenza | stop prima di start, prelievo durante quale linea, date impossibili → warning | S |
| 7. Propagazione | update/extract/remove come per le date dei prelievi | S |

Totale circa **M**.

### Decisioni da prendere
1. **Formato di ingresso:** file separato `treatments.tsv` in formato lungo (una riga per evento) oppure colonne in
   patient.tsv (una per linea)? Consiglio: file separato, più flessibile.
2. **Giorno 0:** la diagnosi (convenzione cBioPortal; serve DIAGNOSIS_DATE in patient.tsv) oppure il primo prelievo
   (come adesso)? Se diagnosi: chi non ce l'ha resta sul primo prelievo, ma nello stesso studio si mescolano due
   riferimenti.
3. **Eventi da supportare ora:** solo TREATMENT, oppure anche STATUS (progressione, risposta, decesso)? E chirurgia /
   radioterapia?
4. **Farmaci:** testo libero o classe di terapia da una lista fissa (la stessa colonna "therapy" di
   resistance_mutations.tsv)? Senza una classe controllata il collegamento con la resistenza (pezzo 5) non funziona.
5. **Risposta RECIST:** la riportiamo solo come evento (valore libero CR/PR/SD/PD) o la validiamo?
6. **OS/PFS (#11):** calcolarli dalla timeline ora o dopo?
7. **Coerenza:** errori bloccanti o solo warning?

### Decisioni #12 (implementata, commit 870c998)
- [x] 1. File separato `treatments.tsv`, formato lungo. Trovato in `[Timeline] TREATMENTS`, altrimenti accanto a
  patient.tsv, altrimenti accanto a sample.tsv (funziona anche con i vendor, es. Guardant).
  Colonne: PATIENT_ID, EVENT_TYPE (TREATMENT/STATUS), START_DATE, STOP_DATE, AGENT, THERAPY_CLASS, LINE, STATUS, RECIST.
- [x] 2. Giorno 0: diagnosi (`DIAGNOSIS_DATE` in patient.tsv, tolta da data_clinical_patient.txt e salvata in
  service/) se ce l'hanno tutti i pazienti con eventi; altrimenti il primo evento di ogni paziente, per tutti, con
  warning. DAYS_FROM_FIRST_SAMPLE resta dal primo prelievo.
- [x] 3. Eventi: TREATMENT e STATUS (Progression, Response, Death; la diagnosi compare come STATUS "Diagnosis").
  Niente chirurgia/radioterapia per ora.
- [x] 4. AGENT testo libero + THERAPY_CLASS da `therapy_classes.tsv` (26 classi, ";" per combinazioni). La colonna
  therapy di resistance_mutations.tsv ora usa le stesse classi (controllato al caricamento).
- [x] 5. RECIST validato: CR/PR/SD/PD, obbligatorio per Response.
- [x] 6. OS/PFS: esclusi.
- [x] 7. Bloccanti: date illeggibili, stop prima di start, campi obbligatori mancanti, valori fuori lista
  (EVENT_TYPE, STATUS, RECIST, THERAPY_CLASS, LINE non numerica). Warning: pazienti non in sample.tsv (ignorati),
  eventi o prelievi dopo il decesso, eventi prima della diagnosi, terapia senza classe, diagnosi mancante per alcuni.

### Output e test
- data_timeline_treatment.txt / data_timeline_status.txt + meta; sezione "Clinical timeline" nel report; colonne
  Compared with / Treatment / Received nelle tabelle del ctDNA seriale.
- Propagazione update/extract/remove: per paziente vince il batch più recente che ha eventi (così una correzione
  sostituisce invece di duplicare).
- Test: errori bloccanti (RECIST, stop<start, classe sconosciuta); due batch uniti = corsa unica; extract ricalcola
  il basale (P1c vs P1a quando P1b manca); senza diagnosi per tutti si usa il primo evento; validatore cBioPortal OK;
  Snakemake identico a varan.py; regressione completa 0 differenze; test seriale precedente identico salvo la nuova
  colonna RESPONSE_BASELINE.

---

## Discussione #4-#7 — Azionabilità, tumore-agnostici, ESCAT, CIViC (implementate, commit 49ac20d)

### Decisioni
- #4 tutto come consigliato: livelli di sensibilità (1-4) e resistenza (R1-R2) separati, Dx/Px no; "azionabile" =
  livelli 1-2 di default, configurabile (`[Actionability] ONCOKB_LEVELS`); CNV e fusioni incluse; versione OncoKB
  salvata (service/ + report). Decisione 6 (senza -k): CIViC è veloce (circa 5 MB, scaricato in meno di 10 s, matching
  in secondi) e senza token, quindi è attivo di default per tutti, con colonne separate da OncoKB.
- #5 tutto come consigliato: gruppo "agnostico" (approvazione) separato da "pan-tumore"; FDA ed EMA in colonne
  distinte con data di verifica; soglia TMB configurabile (10 di default, con nota sul test); stato Detected /
  Not detected / Not tested (gene fuori pannello, nessun dato); una colonna riassuntiva più la tabella nel report.
- #6: Varan fornisce il meccanismo e il formato (escat_template.tsv); la tabella ESCAT del laboratorio va convertita
  in quel formato. AMP/ASCO/CAP solo dalle assertion CIViC, mai calcolato da Varan.
- #7: CIViC affiancato a OncoKB; assertion + evidence A-B; solo annotazione, mai filtro; nessuna segnalazione di
  discordanze.

### Cosa produce
- OncoKB (con -k): ONCOKB_HIGHEST_LEVEL, ONCOKB_RESISTANCE_LEVEL, ONCOKB_N_ACTIONABLE, ONCOKB_ACTIONABLE; livelli CNV
  salvati in service/oncokb_cna_levels.tsv (data_cna.txt non li può contenere); campioni uniti da batch senza -k
  = NA ed esclusi dai conteggi (detto nel report). Riepilogo di coorte per livello e per tipo tumorale.
- CIViC: download automatico al primo uso (o `python civic.py --download`, e nel Dockerfile in fase di build), data
  del download nel report. Matching: cambi proteici, codoni, delezioni/inserzioni di esone, amplificazioni,
  delezioni, fusioni; esclusi i profili generici ("Mutation", espressione...) e le combinazioni. MAF: CIVIC.SENSITIVITY,
  CIVIC.RESISTANCE, CIVIC.EVIDENCE, CIVIC.ASSERTIONS; campione: CIVIC_HIGHEST_LEVEL, CIVIC_RESISTANCE_LEVEL,
  CIVIC_N_ACTIONABLE, CIVIC_ACTIONABLE. Limite dichiarato: le malattie CIViC non sono mappate su OncoTree, quindi i
  livelli valgono "per qualsiasi tumore" e ogni voce riporta la sua malattia.
- Tumore-agnostici (tumor_agnostic.tsv, verificata il 31/10/2026): agnostici NTRK, RET, BRAF V600E, MSI-H, TMB-H;
  pan-tumore ALK, ROS1, FGFR (fusioni, hotspot FGFR3), ERBB2 (amplificazione, mutazioni attivanti), HRD (GIS >= 42).
  HER2 IHC 3+ (T-DXd) non valutabile dal DNA, scritto in nota. Colonne TUMOR_AGNOSTIC_BIOMARKERS,
  TUMOR_AGNOSTIC_NOT_TESTED, PAN_TUMOR_BIOMARKERS.
- ESCAT (`[ESCAT] FILE`): righe gene / alterazione / codici OncoTree (valgono anche per i sottotipi: NSCLC copre LUAD,
  grazie a oncotree_parents.tsv) / livello; tabella sbagliata = la corsa si ferma subito. ESCAT.LEVEL nel MAF,
  ESCAT_HIGHEST_LEVEL e ESCAT_ALTERATIONS per campione.
- Case list: OncoKB azionabili, CIViC A-B, biomarcatore agnostico, MSI-H, TMB-high, Tissue, Plasma.

### Test
- Regressione completa: solo aggiunte (nuove colonne, namespace CIVIC, case list, riepiloghi), nessun dato esistente
  cambiato; validatore cBioPortal OK su create, update, extract, remove, Guardant; Snakemake OK.
- Studio di prova con mutazioni reali e OncoKB simulato: BRAF V600E, EML4::ALK, EGFR esone 19, KRAS G12C/G12D,
  MSI-H e TMB-H riconosciuti; ESCAT applicato per sottotipo (NSCLC → LUAD) e non applicato fuori tumore; tabella
  ESCAT con livello errato blocca la corsa; validatore OK.

### Da fare quando c'è la tabella ESCAT vera
Convertirla nel formato di escat_template.tsv (GENE, ALTERATION, TUMOR, LEVEL, NOTE, SOURCE) e indicarla in
`[ESCAT] FILE`.
