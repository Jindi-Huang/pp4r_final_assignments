import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from stargazer.stargazer import Stargazer
import seaborn as sns

# Histogram of Gross Rent and Gross Rent per Sqm
def histogram_gross_rent(df, output_path):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    ax1.hist(df['Gross Rent (CHF)'], bins=100, alpha=0.5, label='Gross Rent')
    ax1.set_xlabel('Gross Rent (CHF)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Gross Rent')
    ax1.axvline(df['Gross Rent (CHF)'].mean(), color='#00008B', linestyle='dashed', linewidth=1, label='Mean')
    ax1.axvline(df['Gross Rent (CHF)'].median(), color='#FF6347', linestyle='dashed', linewidth=1, label='Median')
    ax2.hist(df['Gross Rent Per Sqm (CHF/m2)'], bins=100, alpha=0.5, label='Gross Rent Per Sqm')
    ax2.set_xlabel('Gross Rent Per Sqm (CHF/m2)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Gross Rent Per Sqm')
    ax2.axvline(df['Gross Rent Per Sqm (CHF/m2)'].mean(), color='#00008B', linestyle='dashed', linewidth=1, label='Mean')
    ax2.axvline(df['Gross Rent Per Sqm (CHF/m2)'].median(), color='#FF6347', linestyle='dashed', linewidth=1, label='Median')
    ax1.legend()
    ax2.legend()
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    plt.close()

# Plot Mean Gross Rent by City, in ascending order
def histogram_gross_rent_by_city(df, output_path):
    df.groupby('City')['Gross Rent (CHF)'].mean().sort_values().plot(kind='bar', alpha=0.5)
    plt.xlabel('City')
    plt.ylabel('Mean Gross Rent (CHF)')
    plt.title('Mean Gross Rent by City')
    plt.savefig(output_path, bbox_inches = "tight")
    plt.close()

# Plot the Distribution of House Characteristics
def histogram_4_characteristics(df, output_path):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
    # Plot the distribution of Living Area
    ax1.hist(df['Living Area (sq. m.)'], bins=100, alpha=0.5, label='Living Area')
    ax1.set_xlabel('Living Area (sq. m.)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Living Area')
    # Plot the distribution of Year Built
    ax2.hist(df['Year Built'], bins=100, alpha=0.5, label='Year Built')
    ax2.set_xlabel('Year Built')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Year Built')
    # Plot the distribution of Number of Rooms, with label step size 0.5
    ax3.hist(df['Number of Rooms'], bins=100, alpha=0.5, label='Number of Rooms')
    ax3.set_xlabel('Number of Rooms')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Number of Rooms')
    ax3.xaxis.set_ticks(np.arange(0, 10, 1))
    # Plot the distribution of Floor Level
    ax4.hist(df['Floor Level'], bins=100, alpha=0.5, label='Floor Level')
    ax4.set_xlabel('Floor Level')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution of Floor Level')
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    plt.close()

# Create a Scatter plot of Gross Rent and Living Area, make the observations in Zurich red, while others blue
def scatter_rent_living_area(df, output_path):
    plt.scatter(df[df['City'] != 'Zürich']['Living Area (sq. m.)'], df[df['City'] != 'Zürich']['Gross Rent (CHF)'], alpha=0.2, color='#00008B')
    plt.scatter(df[df['City'] == 'Zürich']['Living Area (sq. m.)'], df[df['City'] == 'Zürich']['Gross Rent (CHF)'], alpha=0.2, color='#FF6347')
    plt.legend(['Other Cities', 'Zurich'])
    plt.xlabel('Living Area (sq. m.)')
    plt.ylabel('Gross Rent (CHF)')
    plt.title('Gross Rent and Living Area')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()



# Create a scatter plot of Year built and Living area, Year Built and Gross Rent, Year Built and Number of Rooms, Year Built as x-axis, and other variables as y-axis, make the observations in Zurich red, while others blue. The three plots in 3 subplots, with 3 rows and 1 column. Add the mean value of each year. 
def scatter_characteristics_year(df, output_path):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10))
    # Plot the distribution of Living Area
    ax1.scatter(df[df['City'] != 'Zürich']['Year Built'], df[df['City'] != 'Zürich']['Living Area (sq. m.)'], alpha=0.2, color='#00008B')
    ax1.scatter(df[df['City'] == 'Zürich']['Year Built'], df[df['City'] == 'Zürich']['Living Area (sq. m.)'], alpha=0.2, color='#FF6347')
    ax1.set_xlabel('Year Built')
    ax1.set_ylabel('Living Area (sq. m.)')
    ax1.set_title('Living Area and Year Built')
    ax1.axhline(df[df['City'] != 'Zürich']['Living Area (sq. m.)'].mean(), color='#00008B', linestyle='dashed', linewidth=1, label='Mean')
    ax1.axhline(df[df['City'] == 'Zürich']['Living Area (sq. m.)'].mean(), color='#FF6347', linestyle='dashed', linewidth=1, label='Mean')
    ax1.legend(['Other Cities', 'Zurich'])
    # Plot the distribution of Gross Rent
    ax2.scatter(df[df['City'] != 'Zürich']['Year Built'], df[df['City'] != 'Zürich']['Gross Rent (CHF)'], alpha=0.2, color='#00008B')
    ax2.scatter(df[df['City'] == 'Zürich']['Year Built'], df[df['City'] == 'Zürich']['Gross Rent (CHF)'], alpha=0.2, color='#FF6347')
    ax2.set_xlabel('Year Built')
    ax2.set_ylabel('Gross Rent (CHF)')
    ax2.set_title('Gross Rent and Year Built')
    ax2.axhline(df[df['City'] != 'Zürich']['Gross Rent (CHF)'].mean(), color='#00008B', linestyle='dashed', linewidth=1, label='Mean')
    ax2.axhline(df[df['City'] == 'Zürich']['Gross Rent (CHF)'].mean(), color='#FF6347', linestyle='dashed', linewidth=1, label='Mean')
    ax2.legend(['Other Cities', 'Zurich'])
    # Plot the distribution of Number of Rooms
    ax3.scatter(df[df['City'] != 'Zürich']['Year Built'], df[df['City'] != 'Zürich']['Number of Rooms'], alpha=0.2, color='#00008B')
    ax3.scatter(df[df['City'] == 'Zürich']['Year Built'], df[df['City'] == 'Zürich']['Number of Rooms'], alpha=0.2, color='#FF6347')
    ax3.set_xlabel('Year Built')
    ax3.set_ylabel('Number of Rooms')
    ax3.set_title('Number of Rooms and Year Built')
    ax3.axhline(df[df['City'] != 'Zürich']['Number of Rooms'].mean(), color='#00008B', linestyle='dashed', linewidth=1, label='Mean')
    ax3.axhline(df[df['City'] == 'Zürich']['Number of Rooms'].mean(), color='#FF6347', linestyle='dashed', linewidth=1, label='Mean')
    ax3.legend(['Other Cities', 'Zurich'])
    # control the layout of the subplots
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    plt.close()


# Create a correlation matrix of gross rent, additional expenses, living area, number of rooms, floor level, year built
def heat_map(df, output_path):
    df[['Gross Rent (CHF)', 'Additional Expenses (CHF)', 'Living Area (sq. m.)', 'Number of Rooms', 'Floor Level', 'Year Built']].corr()
    # plot the correlation matrix as a heatmap
    
    sns.heatmap(df[['Gross Rent (CHF)', 'Additional Expenses (CHF)', 'Living Area (sq. m.)', 'Number of Rooms', 'Floor Level', 'Year Built']].corr(), annot=True, cmap='coolwarm')
    plt.title('Correlation Matrix')
    # Save heatmap
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

# Regression: What decides Gross Rents?
def regression_table(df, output_path):
    # Preparation: add city dummies, which should be a 0/1 variable for each city, with Zurich as the base city
    df = pd.get_dummies(df, columns=['City'], drop_first=True, dummy_na=True)
    # Select columns starting with "City_"
    city_columns = df.columns[df.columns.str.startswith('City_')]
    # Convert the selected columns to 1 and 0
    df[city_columns] = df[city_columns].astype(int)
    
    model1 = sm.OLS(df['Gross Rent (CHF)'], sm.add_constant(df[['Living Area (sq. m.)']]), missing='drop').fit()
    model2 = sm.OLS(df['Gross Rent (CHF)'], sm.add_constant(df[['Living Area (sq. m.)', 'Number of Rooms']]), missing='drop').fit()
    model3 = sm.OLS(df['Gross Rent (CHF)'], sm.add_constant(df[['Living Area (sq. m.)', 'Number of Rooms', 'Floor Level']]), missing='drop').fit()
    model4 = sm.OLS(df['Gross Rent (CHF)'], sm.add_constant(df[['Living Area (sq. m.)', 'Number of Rooms', 'Floor Level', 'Year Built']]), missing='drop').fit()
    model5 = sm.OLS(df['Gross Rent (CHF)'], sm.add_constant(df[['Living Area (sq. m.)', 'Number of Rooms', 'Floor Level', 'Year Built'] + list(df.columns[df.columns.str.startswith('City_')])]), missing='drop').fit()

    stargazer = Stargazer([model1, model2, model3, model4, model5])
    stargazer.title('Regression on Gross Rent')
    stargazer.covariate_order(['Living Area (sq. m.)', 'Number of Rooms', 'Floor Level', 'Year Built'])
    stargazer.dependent_variable_name('Dependent variable: Gross Rent (CHF)')
    # add a row below to indicate whether we include city dummies
    stargazer.add_line('City FE', ['No', 'No', 'No', 'No', 'Yes'])
    stargazer.show_degrees_of_freedom(False)
    # save table
    html = stargazer.render_html()
    func = open(output_path,"w")
    func.write(html)





def data_analysis(input_path, output_files):
    """Analyze the data and export the results

    Args:
        input_path (str): Path to the processed data (Homegate_data_cleaned.csv)
        output_path (str): Path to the output directory
    """

    df = pd.read_csv(input_path)
    histogram_gross_rent(df, output_files[0])
    histogram_gross_rent_by_city(df, output_files[1])
    histogram_4_characteristics(df, output_files[2])
    scatter_rent_living_area(df, output_files[3])
    scatter_characteristics_year(df, output_files[4])
    heat_map(df, output_files[5])
    regression_table(df, output_files[6])
    


if __name__ == "__main__":   
    data_analysis(
        snakemake.input[0], 
        snakemake.output
    )