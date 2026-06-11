# Obesity Subtype Framework

This repository provides non-core analytical code templates associated with a precision-obesity subtype research framework. The framework is designed to support downstream analyses after obesity subtype labels have already been generated, including disease association analyses, molecular signature construction, and candidate mediation analyses.

This repository contains analytical code only and does not include any individual-level UK Biobank data, participant identifiers, omics data, genetic data, or other restricted participant-level records.

## Repository Scope

The materials in this repository are intended for transparent reporting of downstream statistical workflows. They do not disclose or implement the proprietary obesity subtype assignment framework.

This repository does not include:

- the core obesity subtype assignment algorithm;
- subtype classifier implementation;
- model coefficients, subtype centroids, or scaling parameters;
- online calculator source code;
- UK Biobank individual-level records;
- omics, genetic, imaging, or clinical participant-level datasets;
- any restricted or identifiable participant information.

## Patent Status

The obesity subtype framework, including the core subtype assignment algorithm and related proprietary implementation details, is currently under patent application. Code and materials directly related to subtype generation will be made available after the patent application is successfully completed, subject to institutional, legal, and publication requirements.

## Online Subtype Calculator

An online obesity subtype calculator is available for academic and non-commercial use:

<https://yilinliao520-eng.github.io/obesity-calculator/>

The web tool supports two modes:

- individual calculation using routine clinical measurements;
- batch CSV prediction using `ID`, `SBP`, `HbA1c`, `TG`, `HDL-C`, `Creatinine`, `ALT`, and `WHR`.

The batch module generates a two-column result CSV containing `ID` and `subtype`. According to the web tool description, calculations run locally in the browser and input data are not uploaded to a server.

## Contents

- `scripts/disease_association_cox.R`: Cox proportional hazards model template for subtype-outcome association analyses.
- `scripts/molecular_signature_lasso.R`: LASSO-based molecular signature construction template for a user-specified subtype contrast.
- `scripts/single_feature_mediation.R`: single-feature mediation analysis template for candidate molecular mediators.

## Input Data

The scripts use generic file paths and placeholder variable names. Users should adapt them to their own approved local data structures. Do not commit individual-level data, restricted cohort files, subtype assignment internals, model coefficients, subtype centroids, scaling parameters, or any confidential research materials to this repository.

## License

Academic non-commercial use only. See `LICENSE`.

## Medical Use

This software is provided for academic research use only. It is not intended for clinical diagnosis, treatment decisions, individualized medical advice, or use as a medical device.
