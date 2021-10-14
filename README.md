# COVID-19 Wastewater Surveillance: An Epidemiological Model

## Wastewater-Based Epidemiological Model: A Simple SEIR

Using wastewater surveillance as a continuous pooled sampling technique has been in place in many countries since the early stages of the outbreak of COVID-19. Since the beginning of the outbreak, many research works have emerged, studying different aspects of *viral SARS-CoV-2 DNA concentrations* (viral load) in wastewater and its potential as an early warning method. However, one of the questions that has remained unanswered is the quantitative relation between viral load and clinical indicators such as daily cases, deaths, and hospitalizations. Few studies have tried to couple viral load data with an epidemiological model to relate the number of infections in the community to the viral burden. We propose **a stochastic wastewater-based SEIR model** to showcase the importance of viral load in the early detection and prediction of an outbreak in a community. We built three models based on whether or not they use the case count and viral load data and compared their *simulations* and *forecasting* quality.

We consider a stochastic wastewater-based epidemiological model with four compartments (hidden states) of susceptible (S), exposed (E), infectious (I), and recovered/removed (R) (Fig. \ref{fig:SEIR-scheme}).

.. image:: figs/SEIR_stochastic_col.PNG
