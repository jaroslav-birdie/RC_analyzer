# Verze RC_analyzer

Každá verze je samostatný soubor, aby bylo možné se kdykoli vrátit.
Pojmenování: `RC_analyzer_DDMMYYYY.py` (datum verze); při více verzích za den `_a`, `_b`, …
Vždy pracujeme na **nejnovějším** souboru; starší zůstávají beze změny jako záloha.

| Soubor | Popis změn |
|--------|------------|
| `RC_analyzer_29082025.py` | Původní vydaná verze (29.08.2025). Počítá pouze RC koeficienty (max, A, peak). |
| `RC_analyzer_30062026_a.py` | + Výpočet CRC koeficientů podle Eq.2 (Sunderland et al., J Nucl Med): hodnoty u ROI, vykreslení a export do CSV. RC i CRC se zobrazují současně. |
| `RC_analyzer_30062026_b.py` | + Přepínač zobrazení **RC / CRC** (řeší příliš široké okno). Zobrazuje se vždy jen jedna skupina; popisky, grafy i meze se mění podle režimu, meze se pamatují zvlášť pro RC a CRC. |
| `RC_analyzer_30062026_c.py` | Oprava **exportu** do CSV: pořadí sloupců RC_max, RC_A50, RC_peak, CRC_max, CRC_A50, CRC_peak, diameter, bg_real, bg_measured, bg_diff(%), bg_COV(%); hodnoty pozadí zarovnány do správných sloupců na samostatném řádku bez RC_/CRC_ hodnot. |
| `RC_analyzer_30062026_d.py` | **Rozložení GUI**: MIP obrazy v sekci „MIPs" uspořádány pod sebe (levý nahoře, „up" dole); sekce „Image threshold" vycentrována pod sloupec MIPů. |
| `RC_analyzer_30062026_e.py` | **Rotující MIP** (varianta B): dvojici MIP nahrazuje jediný přední pohled, který se otáčí kolem dlouhé osy fantomu. Sada snímků (36 po 10°, in-plane downsampling faktor 2, `order=0`) se předpočítá hned po segmentaci každého zdroje a rotace se rovnou rozjede; přehrávání jen cyklí hotové obrázky přes `after()` (~10 fps), takže je výkonově nezávislé na stroji. Tlačítko **Spin/Stop** pod obrazem, zoom (Shift+kolečko) zachován. |
| `RC_analyzer_30062026_f.py` | **Předpočet rotace na pozadí** (vlákno) – segmentace zdroje už nečeká na výpočet snímků; během výpočtu se dál točí předchozí výsledek a po dokončení se sada bezešvě prohodí. Generační token zajistí, že při rychlém sledu segmentací vyhraje nejnovější výpočet. Méně projekcí: **24 (po 15°)**; downsampling zvýšen na **faktor 4** kvůli velkým datům (440×440×240 → ~2 s na zdroj). |
| `RC_analyzer_30062026_g.py` | **Adaptivní downsampling** podle velikosti matice (rows×cols): <200 → bez downsamplingu (1), 200–399 → faktor 2, ≥400 → faktor 4. Volí se automaticky při načtení dat. |
| `RC_analyzer_30062026_h.py` | **Reset zoomu MIP** prostředním tlačítkem myši (stejně jako u zobrazení fantomu): vrátí zvětšení i posun rotujícího MIPu do výchozího stavu. **Aktuální pracovní verze.** |

> Pozn.: kompletní historie změn je i v gitu (`git log`), tato tabulka slouží pro rychlou orientaci v souborech.
