import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import resample
import numpy as np

# Read the data from a TSV file
df = pd.read_csv('/data/san/data0/users/david/intelligence/tables/ctrl_v_LanM_vs_Pks.tsv', sep='\t')

# Print the columns to check for the presence of 'pfam'
print("Columns in the DataFrame:", df.columns)

# Convert 'NA' to None for proper handling
df['pfam'] = df['pfam'].replace('NA', None)

# Create a pivot table
pivot_df = df.pivot_table(index='Nucleotide_acc', columns='pfam', aggfunc='size', fill_value=0)

# Print head of pivot
print("Head of pivot table:")
print(pivot_df.head())

# Reset index to bring 'Nucleotide_acc' back as a column
pivot_df.reset_index(inplace=True)

# Merge with the original dataframe to get the 'group' labels back
df_group = df[['Nucleotide_acc', 'group']].drop_duplicates()
merged_df = pd.merge(pivot_df, df_group, on='Nucleotide_acc')

# Encoding the 'group' column
le = LabelEncoder()
merged_df['group'] = le.fit_transform(merged_df['group'])

# Define features and labels
X = merged_df.drop(columns=['Nucleotide_acc', 'group'])
y = merged_df['group']

# Train a Random Forest classifier to determine feature importances
clf_initial = RandomForestClassifier(n_estimators=100, random_state=42)
clf_initial.fit(X, y)

# Determine feature importance
initial_feature_importances = clf_initial.feature_importances_
initial_feature_names = X.columns

# Create a DataFrame for the feature importances
initial_importance_df = pd.DataFrame({
    'Feature': initial_feature_names,
    'Importance': initial_feature_importances
})

# Sort the DataFrame by importance
initial_importance_df = initial_importance_df.sort_values(by='Importance', ascending=False)

# Identify the top 50 most important features
top_50_features = initial_importance_df.head(50)['Feature'].tolist()

# Remove the top 50 most important features from the dataset
X_filtered = X.drop(columns=top_50_features)

# Split the filtered data into training, validation, and testing sets
X_train_filtered, X_temp_filtered, y_train, y_temp = train_test_split(X_filtered, y, test_size=0.3, random_state=42)
X_val_filtered, X_test_filtered, y_val, y_test = train_test_split(X_temp_filtered, y_temp, test_size=0.5, random_state=42)

# Train a Random Forest classifier on the filtered dataset
clf_filtered = RandomForestClassifier(n_estimators=100, random_state=42)
clf_filtered.fit(X_train_filtered, y_train)

# Make predictions on the test data
y_pred_filtered_test = clf_filtered.predict(X_test_filtered)
y_pred_filtered_val = clf_filtered.predict(X_val_filtered)

# Evaluate the model
accuracy_filtered_test = accuracy_score(y_test, y_pred_filtered_test)
accuracy_filtered_val = accuracy_score(y_val, y_pred_filtered_val)
report_filtered_test = classification_report(y_test, y_pred_filtered_test)
report_filtered_val = classification_report(y_val, y_pred_filtered_val)

print(f"Accuracy on test set after filtering: {accuracy_filtered_test}")
print("Classification Report on test set after filtering:")
print(report_filtered_test)

print(f"Accuracy on validation set after filtering: {accuracy_filtered_val}")
print("Classification Report on validation set after filtering:")
print(report_filtered_val)

# Determine feature importance on the filtered dataset
filtered_feature_importances = clf_filtered.feature_importances_
filtered_feature_names = X_filtered.columns

# Create a DataFrame for the feature importances on the filtered dataset
filtered_importance_df = pd.DataFrame({
    'Feature': filtered_feature_names,
    'Importance': filtered_feature_importances
})

# Sort the DataFrame by importance
filtered_importance_df = filtered_importance_df.sort_values(by='Importance', ascending=False)

# Save the feature importance table
filtered_importance_df.to_csv('filtered_feature_importance.csv', index=False)

# Print the most important features after filtering
print("Top 10 most important features after filtering:")
print(filtered_importance_df.head(10))

# Confusion matrix
cm = confusion_matrix(y_test, y_pred_filtered_test)
plt.figure(figsize=(10, 7))
sns.heatmap(cm, annot=True, fmt='g')
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix - Test Set')
plt.savefig('confusion_matrix_test.png')
plt.show()

# Bootstrap replicates for the whole set
n_iterations = 1000
n_size = int(len(X_filtered) * 0.7)
bootstrap_accuracies = []

for _ in range(n_iterations):
    # Create a bootstrap sample
    X_sample, y_sample = resample(X_filtered, y, n_samples=n_size, random_state=42)
    # Train the model
    clf_bootstrap = RandomForestClassifier(n_estimators=100, random_state=42)
    clf_bootstrap.fit(X_sample, y_sample)
    # Evaluate the model on the test set
    y_pred_bootstrap = clf_bootstrap.predict(X_test_filtered)
    accuracy = accuracy_score(y_test, y_pred_bootstrap)
    bootstrap_accuracies.append(accuracy)

# Plot accuracy distribution
plt.figure(figsize=(10, 7))
plt.hist(bootstrap_accuracies, bins=30, edgecolor='k', alpha=0.7)
plt.axvline(np.mean(bootstrap_accuracies), color='r', linestyle='dashed', linewidth=1)
plt.title('Bootstrap Accuracy Distribution')
plt.xlabel('Accuracy')
plt.ylabel('Frequency')
plt.savefig('bootstrap_accuracy_distribution.png')
plt.show()
