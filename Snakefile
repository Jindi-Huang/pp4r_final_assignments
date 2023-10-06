import platform


# Create website
rule create_website:
    input:
        "output/figures/histogram_gross_rent.png",
        "output/figures/histogram_gross_rent_by_city.png",
        "output/figures/histogram_4_characteristics.png",
        "output/figures/scatter_rent_living_area.png",
        "output/figures/scatter_characteristics_year.png",
        "output/figures/heat_map.png",
	"output/figures/regression_table.html"
    output:
        "website.html"

# Data analysis
rule data_analysis:
    conda: "envs/scraping.yaml"
    input:
        "data/Homegate_data_cleaned.csv"
    output:
        "output/figures/histogram_gross_rent.png",
	"output/figures/histogram_gross_rent_by_city.png",
        "output/figures/histogram_4_characteristics.png",
        "output/figures/scatter_rent_living_area.png",
        "output/figures/scatter_characteristics_year.png",
        "output/figures/heat_map.png",
	"output/figures/regression_table.html"
    script:
        "src/homegate_analysis.py"



# Data cleaning
rule data_cleaning:
    conda: "envs/scraping.yaml"
    input:
        "data/Homegate_data.csv"
    output:
        "data/Homegate_data_cleaned.csv"
    script:
        "src/homegate_cleaning.py"



# Data scraping
rule data_scraping:
    conda: "envs/scraping.yaml"
    input:

    output:
        "data/Homegate_data.csv"
    script:
        "src/homegate_scraping.py"
