import platform

# Create website
rule create_website:
    input:
        "output/figures/histogram_gross_rent.png"
        "output/figures/histogram_gross_rent_by_city.png"
        "output/figures/histogram_4_characteristics.png"
        "output/figures/scatter_rent_living_area.png"
        "output/figures/scatter_characteristics_year.png"
        "output/figures/heat_map.png"
    output:
        "website.html"

# Postcode analysis
rule postcode_analysis:
    input:
        "data/Homegate_data_cleaned.csv"
    output:
        "data/"
    script:
        "src/homegate_postalcodes.ipynb"

# Data analysis
rule data_analysis:
    input:
        "data/Homegate_data_cleaned.csv"
    output:
        "output/figures/histogram_gross_rent.png"
        "output/figures/histogram_gross_rent_by_city.png"
        "output/figures/histogram_4_characteristics.png"
        "output/figures/scatter_rent_living_area.png"
        "output/figures/scatter_characteristics_year.png"
        "output/figures/heat_map.png"
    script:
        "src/homegate_analysis.ipynb"

# Data cleaning
rule data_cleaning:
    input:
        "data/Homegate_data.csv"
    output:
        "data/Homegate_data_cleaned.csv"
    script:
        "src/homegate_cleaning.py"
