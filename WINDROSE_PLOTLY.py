import plotly.express as px
df = px.data.wind()
print(df)
df.to_csv("w.csv")
fig = px.bar_polar(df, r="frequency", theta="direction",
                   color="strength", template="plotly_dark",
                   color_discrete_sequence= px.colors.sequential.Plasma_r)
fig.show()