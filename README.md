# Obesity Subtype Framework

This repository contains non-core analysis templates for the obesity subtype research project. The scripts are intended to support reproducible downstream analyses after subtype assignments have already been generated.

The repository does not include the core obesity subtype assignment algorithm, subtype classifier, model parameters, subtype centroids, scaling parameters, or online calculator source code.

## Contents

- `scripts/disease_association_cox.R`: Cox proportional hazards models for subtype-outcome associations.
- `scripts/molecular_signature_lasso.R`: LASSO-based molecular signature construction for a user-specified subtype contrast.
- `scripts/single_feature_mediation.R`: single-feature mediation analysis for candidate molecular mediators.

## Input Data

The scripts use generic placeholders and should be adapted to local data structures. Individual-level data, restricted cohort files, subtype assignment internals, model coefficients, centroids, and scaling parameters should not be committed to this repository.

## License

Academic non-commercial use only. See `LICENSE`.

## Medical Use

This software is provided for academic research use only and is not intended for clinical diagnosis, treatment decisions, individualized medical advice, or use as a medical device.
