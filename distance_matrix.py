from pipeline import *

if __name__ == "__main__":
    getDissimilaritiesMatrix('T min à 2m (C)', "output/meteo_data_matrix/t_min_matrix")
    getDissimilaritiesMatrix(
        'T max à 2m (C)', "output/meteo_data_matrix/t_max_matrix")
    getDissimilaritiesMatrix('Humidité relative à 2m(%)',
                             "output/meteo_data_matrix/humidity_matrix")
    getDissimilaritiesMatrix('Précipitation totale sur le mois (mm)',
                             'output/meteo_data_matrix/precipitation_matrix')
    getDissimilaritiesMatrix('Pression en surface (kPa)',
                             'output/meteo_data_matrix/pressure_matrix')

