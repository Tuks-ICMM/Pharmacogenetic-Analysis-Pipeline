__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Sarah Turner",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"

# %%

import plotly.express as px

# %%
df = px.data.iris().rename(
    columns={
        "sepal_width": "Sepal Width",
        "sepal_length": "Sepal Length",
        "petal_length": "Petal Length",
        "species": "Species",
        "petal_width": "Petal Width",
    }
)  # iris is a pandas DataFrame
fig = px.scatter(
    df,
    x="Sepal Width",
    y="Sepal Length",
    color="Species",
    symbol="Species",
    hover_data=["Petal Width"],
    marginal_x="box",
    marginal_y="box",
    title="Population Stratification PCA",
    color_discrete_sequence=px.colors.qualitative.Bold
)
# fig.update_traces(marker_size=10)
fig.update_layout(
    scattermode="group"
)  # Set the overlay when two points co-occur at same coordinates
fig.show()

# %%
