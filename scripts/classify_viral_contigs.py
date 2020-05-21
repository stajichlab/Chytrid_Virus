#!/usr/bin/env python3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import model_selection

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report, confusion_matrix

from sklearn.model_selection import RandomizedSearchCV

dfall = pd.read_csv("results/genome_stats_nofilter/Burma_1F.sum_stat.tsv", delimiter="\t")
print(dfall.shape)

df = dfall.drop(['ID','LENGTH'],axis=1)
print(df.shape)
X = df.drop('CAPSID', axis=1)
y = df['CAPSID']
print(y.shape)
# implementing train-test-split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=66)

print(X_train,X_test,y_train,y_test)
# random forest model creation
rfc = RandomForestClassifier(n_estimators=100, warm_start=True)
rfc.fit(X_train,y_train)
# predictions
rfc_predict = rfc.predict(X_test)

rfc_cv_score = cross_val_score(rfc, X, y, cv=10, scoring='roc_auc')

print("=== Confusion Matrix ===")
print(confusion_matrix(y_test, rfc_predict))
print("\n")
print("=== Classification Report ===")
print(classification_report(y_test, rfc_predict))
print("\n")
print("=== All AUC Scores ===")
print(rfc_cv_score)
print("\n")
print("=== Mean AUC Score ===")
print("Mean AUC Score - Random Forest: ", rfc_cv_score.mean())

# number of trees in random forest
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# number of features at every split
max_features = ['auto', 'sqrt']

# max depth
max_depth = [int(x) for x in np.linspace(100, 500, num = 11)]
max_depth.append(None)
# create random grid
random_grid = {
 'n_estimators': n_estimators,
 'max_features': max_features,
 'max_depth': max_depth
 }
# Random search of parameters
rfc_random = RandomizedSearchCV(estimator = rfc, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
# Fit the model
rfc_random.fit(X_train, y_train)
# print results
print(rfc_random.best_params_)