## Background

The original R script performs the following steps:

1. **Data Import & Cleaning**  
2. **Risk Model Training & Validation**  
3. **Score Computation**  
4. **Export of Results** (CSV, JSON)

This pipeline is used to generate patient‐level risk scores based on clinical and genomic covariates. By converting to Python, we can:

- Use **pandas** & **scikit-learn** (or **statsmodels**) for data handling and modeling.
- Integrate with modern ML frameworks (e.g., **XGBoost**, **LightGBM**).
- Deploy on cloud platforms using Docker, Kubernetes, or serverless functions.
