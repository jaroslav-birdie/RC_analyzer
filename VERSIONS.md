# Verze RC_analyzer

Každá verze je samostatný soubor, aby bylo možné se kdykoli vrátit.
Pojmenování: `RC_analyzer_DDMMYYYY.py` (datum verze); při více verzích za den `_a`, `_b`, …
Vždy pracujeme na **nejnovějším** souboru; starší zůstávají beze změny jako záloha.

| Soubor | Popis změn |
|--------|------------|
| `RC_analyzer_29082025.py` | Původní vydaná verze (29.08.2025). Počítá pouze RC koeficienty (max, A, peak). |
| `RC_analyzer_30062026_a.py` | + Výpočet CRC koeficientů podle Eq.2 (Sunderland et al., J Nucl Med): hodnoty u ROI, vykreslení a export do CSV. RC i CRC se zobrazují současně. |
| `RC_analyzer_30062026_b.py` | + Přepínač zobrazení **RC / CRC** (řeší příliš široké okno). Zobrazuje se vždy jen jedna skupina; popisky, grafy i meze se mění podle režimu, meze se pamatují zvlášť pro RC a CRC. **Aktuální pracovní verze.** |

> Pozn.: kompletní historie změn je i v gitu (`git log`), tato tabulka slouží pro rychlou orientaci v souborech.
