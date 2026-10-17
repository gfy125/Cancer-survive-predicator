# Skin Cancer Survival Time Calculator

Skin cancer is the abnormal growth of skin cells, often caused by excessive exposure to ultraviolet (UV) radiation from the  sun or tanning beds. It is the most common type of cancer and typically develops on areas of the skin that are frequently exposed to sunlight, such as the face, neck, and arms. There are three main types: basal cell carcinoma, squamous cell carcinoma, and melanoma, the most dangerous form, which can spread to other parts of the body if not detected early. Preventative measures include using sunscreen, wearing protective clothing, and avoiding sun exposure during peak hours.  Early detection and treatment significantly improve outcomes.

This calculator utilises melanoma patient demographics and gene expression data to predict survival probability and provides healthcare professionals with personalised insights that may help improve patient outcomes and overall well-being.

## How To Use the App
All the following will take place in the skin_cancer_risk folder.

### How to Obtain the Data
Run cox_jeslyn_trials_predicting.qmd to obtain the model.

### Running the App
Run the App from app.R

### Predicting
Follow the instructions on the 'About' tab or:
1. Navigate to the Calculator tab.
2. Input patient details, including:
    1. Age
    2. Gender
    3. Cancer Stage
3. Upload the patient's gene expression data (.csv or .RData file) or use default values.
4. Click on the Predict button.
5. View the results including survival curve, risk assessment, and recommended lifestyle changes.

### Other Features
Other than predicing a patients survival time. Users can also view sections of the Exploratory Data Analysis in the 'Model Explanation' tab, this explains the most features creating most risk and protection, the risk score distribution, and the Concordance of the model.

Another feature is the 'Dataset Explorer' where users are given the summary statistics and are able to explore the clinical and gene expression data that was used to train the model.
