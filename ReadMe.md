# biocharStability - dataset of biochar incubations, alongside code to analyse it

|           **Info** 	| Dataset of biochar incubations, alongside code to analyse it	|
|-------------------:	|---------------------------------------------	|
|        **Website** 	| https://biochar.systems/stability/  <br> https://biochar.systems/durability-statement/       	|
|        **Webinars** 	| https://www.youtube.com/watch?v=FqWrpOK8RS0 (2023-04 / EGU)	<br> https://www.youtube.com/watch?v=VXL3c4kX9tY (2022-06 / NegCO2) |
|    **Start / End** 	| 2022-01 to 2025                             	|
|   **Contact points** 	| kolinlagringattraknamed@2050.se (general request) <br> elias@ecoleaf.consulting (data & library) |
|       **Publication** 	| *manuscript under submission, pre-print not yet available*                    	|
|       **How to cite?** 	|  Azzi, E. S. (2023). biocharStability - dataset of biochar incubations, alongside code to analyse it [Dataset & Code]. https://github.com/SLU-biochar/biocharStability                   	|
|       **License** 	| [![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]                             	|

The climate change mitigation benefits of biochar systems arise to a large extent from **carbon storage in biochar** However, there are on-going scientific discussions on how to precisely estimate the **durability of biochar carbon storage in soils**. Estimates vary from decades to millennia, building on different evidence, modelling approaches and theories.

In this project (2022-2025), [funded by the Swedish Energy Agency](https://bioplusportalen.se/en/project/biochar-stability-validation-reaching-a-new-level-of-understanding-and-transparency), the knowledge is strengthened by performing:
- novel analysis of available data from biochar incubation experiments (started 2022-03, ended 2023-07)
- new incubation trials (started 2023-02)
- new field trials (started 2023-05)

The project will involve Swedish and international biochar stakeholders and will develop guidelines for estimation of biochar stability based on biochar origin, properties, and use. Since this topic requires long-term research, field trials are established, planned for long-term analysis, as well as a web app for analysis of future research results.

The project can be important for the development of biochar in voluntary markets and policy for carbon dioxide removal, also known as negative emissions.

*Vocabulary remark: the word ``stability`` is used in the name of the library, as it was still a common term back when the project started; however, these days, words like ``persistence``, ``permanence``, or ``durability`` are preferred depending on which background one has.*

## This library ``biocharStability``

This library contains both the ``dataset`` compiled and the ``code`` developed during the project to analyse biochar incubation experiments, and model durability of carbon storage with biochar.

- The ``database`` is available as a .xlsx file, in `` biocharStability / database / `` alongside other relevant data (e.g. Q<sub>10</sub> values). A description of the database structure is provided in  `` biocharStability / database / schema ``. Simply open the html file in your favorite web browser.
- Several ``former assessments`` that have released data in Excel are compiled under `` biocharStability / database / former-assessment /`` 
- The ``code`` is provided in the folder `` biocharStability `` and is split in a handful of python scripts that reflects the different parts of the modelling, e.g. analyse.py, dashboard.py, visualize.py, utils.py 
- Demonstration and development ``notebooks `` are provided:
    - in the folder `` notebooks-demo ``, curated notebooks are shared to give an introduction of how to work with the library. These notebooks are a good starting point for new users.
    - in the folder `` notebooks-dev ``, notebooks used during development and  research are shared. These notebooks can be interesting for advanced users. Warning: they are provided without guarantee.  
- Notebooks used for the academic publication associated to this work are in the folder ``manuscript`` (not yet available)
- Outputs of simulations can be saved in a folder ``simulation`` that is not synced on GitHub (cf. gitignore)

## HOW TO INSTALL OR USE?

### as a ``python`` user / developper / contributor

1. Fork / Clone the repository to your local machine

2. Create a development environment (e.g. in conda ``` conda create -n biocharstab python=3.9 ```)

3. In this environment, from command line, change directory to the forked/cloned repository, and ``` pip install -e . ```

    Note: "-e" stands for --editable installation

    Note: "." stands for current directory

    Note: the setup.py file of this package contains a list of required packages that will also be installed. If anything is missing or if there are version clashes, please check also the file ``environment.yml`` in the root folder. Ask for helf, if needed.

4. In your python shell, notebook or development environment, load the now pip-installed libray with ``` import biocharStability as bs ```

    Note: example notebooks are also available in this repository, for a quick start.

5. Check the documentation available in the folder ``docs``: simply open the file ``index.html`` in your favorite web browser.

    Memo on how to update the documentation, built with pandoc:
        
        a. From command line, activate conda environment and cd to folder ``../biocharStability``
        b. Use: `` pdoc biocharStability # `` 
        c. Use: `` pdoc biocharStability -o docs `` to save in folder ``/docs`` the html file

### as a ``non-python`` user

You can:

1. Explore the files, here on Github or by downloading the whole repository.
For instance, you can check the ``database`` as an Excel file, which you will find in `` biocharStability / database / * here * .xlsx ``. Other Excel files and documentation is also available.

2. Check out the data and some analyses on the website of the website, more specifically, its data dashboard here https://biochar.systems/stability/app (only partial data overview)

## CONTRIBUTION GUIDE

Contributions in any form are welcome! If you have any doubt, [don't hesitate to e-mail us](mailto:elias@ecoleaf.consulting?subject=[Contribution%20to%20Biochar%20Stability%20Project])

### 1/ DATA CONTRIBUTIONS

You can for instance:

- Provide incubation data and metadata for new observations (in any format .csv, .xslx, .json or similar), whether or not the data is formatted according to our database structure
    - Technical note: We intended *primarily* to collect biochar incubation data of experiments that can distinguish biochar carbon flows from soil/biomass carbon flows (e.g. via isotopic technines). However, we are *also* attempting to collect biochar incubation data where only total carbon flows were measured, provided data from a control setup is also available.

- Provide corrections, additions, validation elements for existing observations
    - Technical note: We intended *primarily* to collect the biochar carbon fluxes, as these were the most ready available. However, we are *also* attempting to collect all timeseries from the experiments (e.g. total carbon flux, soil carbon fluxes, soil carbon fluxes in control case, in absolute and relative units, isotopic fraction measurements, temperature, moisture...). Much help is needed there.
    - Technical note: data validation or verification is also important - anyone can verify the information contained in the ``metadata`` table, and update accordingly the ``validation`` table.

- Recommend new articles for our consideration that are not yet listed in the ``articles`` table.

### 2/ CODE CONTRIBUTIONS

You can for instance:

- Improve the data handling code in ``utils.py`` (mainly relying on numpy and pandas)
- Improve the data visualisation code in ``visualize.py`` (mainly relying on matplotlib, seaborn, ternary)
- Improve the dashboard code in ``dashboard.py`` (mainly relying on bokeh, with exports as .js components integrated to a markdown app.md file for our static website in hugo available at https://biochar.systems/stability/app)
- Improve the analysis code in  ``analyse.py`` (mainly relying on numpy, pandas, scipy.optimize, statsmodels, scikit-learn)
- Suggest & develop new functionalities or analyses in notebooks or custom scripts

### 3/ OTHER CONTRIBUTIONS

You can for instannce:
- Develop or improve the documentation, correct typos
- Suggest new functionalities or analyses
- Share notebooks or work you have done using this repository :blush:

For any contribution, [just email us](mailto:elias@ecoleaf.consulting?subject=[Contribution%20to%20Biochar%20Stability%20Project]). Thank you for your efforts.

## LICENSE

[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
Attribution-ShareAlike 4.0 International


## Documentation built with ``pandoc``

