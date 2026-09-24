# Data

## Kaggle Wind Turbine SCADA (`raw/T1.csv`)

- Source: https://www.kaggle.com/datasets/berkerisen/wind-turbine-scada-dataset (version 3,
  last updated 2019-03-07). Licence listed as "Unknown" on Kaggle: cite, do not redistribute.
- Downloaded 2026-09-24 via `kagglehub.dataset_download("berkerisen/wind-turbine-scada-dataset")`
  (anonymous download works; no API key needed). Raw file is git-ignored.
- SHA-256: `f92c33c1cc199756b759fec251eeddc506f35d6a9e47ede034fc3e339a7c46d5`
- Single turbine in Turkey, 2018-01-01 00:00 to 2018-12-31 23:50, 10-minute resolution.
- Columns: `Date/Time` (format `%d %m %Y %H:%M`, file starts with a UTF-8 BOM),
  `LV ActivePower (kW)`, `Wind Speed (m/s)`, `Theoretical_Power_Curve (KWh)`, `Wind Direction (°)`.

### Facts checked on download

- 50,530 rows; full 10-min grid would be 52,560, so 2,030 timestamps missing, no duplicates.
  Largest gaps: 625 slots from 2018-01-26 06:20, 546 from 2018-09-28 21:20, 520 from 2018-11-10
  21:10; the rest are small (73 slots and fewer).
- Speed 0 to 25.2 m/s; 10 exact zeros, 121 below 0.5 m/s (matters for `log(speed)`).
- Direction 0 to 360 degrees. Two main directions: about 60 degrees (56% of obs within 0-90) and about 200
  degrees (20% within 180-225). Circular mean 65 degrees, mean resultant length 0.29.
- Median absolute 10-min direction change 2.7 degrees; 99th percentile 62 degrees.
- Wraps between consecutive observations: 696 with the (-pi, pi] cut, 709 with [0, 2pi); the fewest
  is 454, with the cut at 130 degrees.
