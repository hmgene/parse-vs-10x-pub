import pandas as pd
from mlxtend.frequent_patterns import apriori, association_rules

# Sample data (replace this with your actual DataFrame)
df = pd.DataFrame({
    'group': ['UH14', 'UH14', 'UH16', 'UH16', 'UH17'],
    'diagnosis': ['GBM', 'GBM', 'GBM', 'GBM', 'Meningioma'],
    'celltype': ['T-cell', 'B-cell', 'Myeloid', 'T-cell', 'B-cell']
})

j= ["diagnosis","leiden"]
df_encoded = pd.get_dummies(tt.obs[j])

frequent_itemsets = apriori(df_encoded, min_support=0.01, use_colnames=True)
rules = association_rules(frequent_itemsets, metric="lift", min_threshold=0.5)

# Display rules sorted by confidence
print(rules[['antecedents', 'consequents', 'support', 'confidence', 'lift']].sort_values(by='confidence', ascending=False))

