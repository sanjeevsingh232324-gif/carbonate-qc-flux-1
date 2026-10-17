Analyze the ocean carbonate chemistry observations in `/app/bottle.csv` and estimate monthly net air–sea CO₂ exchange from `/app/underway.csv`. Save the final result as valid JSON at `/app/report.json`.
The bottle dataset contains `sample_id`, `station`, `lat`, `lon`, `depth_m`, `pressure_dbar`, `temp_c`, `salinity`, `dic_umolkg`, `ta_umolkg`, `ph_meas`, `ph_temp_c`, `fco2_uatm`, `phosphate_umolkg`, and `silicate_umolkg`. Temperatures are in °C, pressure is in dbar, salinity is dimensionless, fCO₂ is in µatm, and DIC, TA, phosphate, and silicate are in µmol/kg. The underway dataset contains `timestamp`, `sst_c` (SST in °C), `sss`, `fco2_sw_uatm` (seawater fCO₂ in µatm), `xco2_atm_ppm` (atmospheric dry-air xCO₂ in ppm), `slp_hpa` (sea-level pressure in hPa), and `wind_10m_ms` (10 m wind speed in m/s). In the bottle data, `ph_meas` is pH on the total scale measured at `ph_temp_c` and 1 atm, and `fco2_uatm` is referenced to the in-situ temperature `temp_c` at 1 atm.
Identify bottle samples with faulty or inconsistent carbonate chemistry observations. Assign each flagged sample exactly one of these diagnosis labels:

* `ph_seawater_scale` — correct `ph_meas`.
* `ph_temp_mislabeled` — correct `ph_temp_c`.
* `dic_ta_swapped` — correct `dic_umolkg` and `ta_umolkg`.
* `dic_volumetric_units` — correct `dic_umolkg`.
* `ta_volumetric_units` — correct `ta_umolkg`.

Use the labels exactly as written. Base diagnoses and corrections on the observed data and their consistency with the other measurements. Apply only corrections supported by the observations, preserve sample IDs, and leave unaffected variables unchanged. For volumetric-unit corrections, convert µmol/L to µmol/kg using seawater density evaluated at 25 °C and 1 atm at the sample's salinity. Do not flag observations without sufficient evidence.
Use Lueker et al. (2000) for the carbonic acid dissociation constants, with pressure corrections from Millero (1995); Dickson (1990) for the bisulfate dissociation constant; Lee et al. (2010) for total boron; and Dickson and Riley (1979) for the fluoride dissociation constant. Report pH on the total scale. Apply consistent carbonate chemistry conventions when evaluating the observations.
Estimate net air–sea CO₂ exchange using the Weiss (1974) CO₂ solubility formulation and the Wanninkhof (2014) gas-transfer parameterization. Convert solubility from mol kg⁻¹ atm⁻¹ to mol m⁻³ atm⁻¹ using seawater density. Calculate atmospheric CO₂ fugacity from dry-air xCO₂ and sea-level pressure, using Weiss and Price (1980) for water-vapour pressure and Weiss (1974) for the fugacity correction, including its virial coefficients. Use the gas-transfer coefficient (k=0.251U_{10}^{2}(Sc/660)^{-0.5}), where (k) is in cm h⁻¹ and (Sc) is the seawater Schmidt number calculated using the Wanninkhof (2014) polynomial. Convert all quantities consistently, including the gas-transfer coefficient, when calculating fluxes. Underway seawater fCO₂ is already referenced to SST.
Define flux as positive for ocean-to-atmosphere transfer (outgassing) and negative for atmosphere-to-ocean transfer (uptake). Calculate the mean flux rate for each calendar month and multiply it by the number of hours in that month. Report monthly integrated fluxes in mol m⁻² month⁻¹.
Write the output JSON with exactly four top-level keys:

* `flagged`: an array of flagged sample IDs.
* `diagnosis`: a mapping from each flagged sample ID to exactly one permitted diagnosis label.
* `corrected`: a mapping from sample IDs to objects containing only the corrected variable names and their numerical values. For example, `{ "sample_id": { "ph_meas": 8.05 } }`.
* `monthly_flux`: a mapping from calendar month to its numerical integrated flux. If the data cover one calendar year, use string keys `"1"` through `"12"`. If multiple years are present, use `"YYYY-MM"` keys.

Ensure the JSON is valid, contains no additional top-level keys, uses the specified field names and units, and contains only finite numerical values. Keep the diagnosis mapping consistent with `flagged` and the corrections supplied. Corrections are checked against reference values with these absolute tolerances: `ph_meas` 0.003, `ph_temp_c` 0.01, `dic_umolkg` and `ta_umolkg` 0.5 µmol/kg for a swap and 1.5 µmol/kg for a volumetric-unit correction. Monthly fluxes are checked to within 1 percent (or 0.0005 mol m⁻² month⁻¹, whichever is larger).
You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
