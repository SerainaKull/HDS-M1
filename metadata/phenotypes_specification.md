# Computational Phenotype Specification: LACDHS Diabetic Retinopathy ARD (v1.0.0)

- **PHENO-01 (`dr_status`):** Binary outcome (`1` = Diabetic Retinopathy present at teleretinal/ophthalmologic exam; `0` = absent).
- **PHENO-02 (`glycemic_control_cat`):** Operational HbA1c groups: <7%; >=7% and <9%; >=9%. Missing cleaned HbA1c produces a missing category. Internal labels do not establish individual treatment targets.
- **PHENO-03 (`stage2_htn_flag`):** TRUE if cleaned SBP >=140 mmHg or DBP >=90 mmHg; NA if either pressure is missing. This is an operational threshold flag, not a hypertension diagnosis.
- **PHENO-04 (`renal_azotemia_flag`):** TRUE if cleaned BUN >20 mg/dL; NA if BUN is missing. This flag does not establish renal disease.
- **PHENO-05 (`miss_*` flags):** Missingness after cleaning, including original missingness and values newly set to NA. lab_missing_count sums six grouped indicators; BP counts once if either pressure is missing. These flags do not establish why data are missing.
