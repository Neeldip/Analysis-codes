import pandas as pd
df = pd.read_excel(r"I:\asos\GHY\ghyapril2018.xlsx",sheet_name="Sheet4")
print(df)
#df = df.to_string
df[(df['time'] =='*:30:00')]
print(df)