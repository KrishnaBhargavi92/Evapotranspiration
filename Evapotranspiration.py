#!/usr/bin/env python3
#Author: P Krishna Bhargavi


import numpy as np

class NetRadiationComponents:

    def __init__(self, **kwargs):       
        if 'albedo' in kwargs: self._albedo = kwargs['albedo']              
        if 'air_temperature_data' in kwargs: self._air_temperature_data = kwargs.get('air_temperature_data')        
        if 'insolation' in kwargs: self._insolation = kwargs.get('insolation')

    def NetShortWaveRadiation(self):
        return(self._insolation*(1-self._albedo))

    def IncomingLongWaveRadiation(self):
        stephan_boltzmann_constant = 5.67*pow(10,-8)
        air_emissivity = 1 #For grey body
        return(stephan_boltzmann_constant*air_emissivity*pow(self._air_temperature_data,4))

class NetRadiation(NetRadiationComponents):

    def __init__(self, **kwargs):
        if 'out_going_long_wave_radiation' in kwargs: self._out_going_long_wave_radiation = kwargs.get('out_going_long_wave_radiation')
        if 'elevation' in kwargs: self._elevation = kwargs.get('elevation')
        super().__init__(**kwargs)

    def NetRadiationValue(self):
        ilr = NetRadiationComponents.IncomingLongWaveRadiation(self)
        nswr = NetRadiationComponents.NetShortWaveRadiation(self)
        return(nswr+(ilr-self._out_going_long_wave_radiation))

class SoilHeatFlux(NetRadiation):

    def __init__(self, **kwargs):
        if 'land_surface_temperature' in kwargs:self._land_surface_temperature = kwargs.get('land_surface_temperature')
        if 'NDVI' in kwargs: self._NDVI = kwargs.get('NDVI')
        super().__init__(**kwargs)

    def SoilHeatFluxValue(self):
        net_radiation_value = NetRadiation.NetRadiationValue(self)
        return(net_radiation_value*((self._land_surface_temperature/self._albedo)*((0.0038*self._albedo)+(0.0074*pow(self._albedo,2)))*(1-0.98*pow(self._NDVI,4)))) ##https://www.sciencedirect.com/science/article/pii/S0895717710005303 --equation 3

class SlopeOfSaturatedVaporPressure():

    def delta(self): #Slope of saturation vapor pressure is represented as delta by Modified Priestley Taylor method
        return((2503.0580/(237.3000+self._air_temperature_data))*np.exp((17.2700*self._air_temperature_data)/(self._air_temperature_data+237.3000)))

class PsychometricConstant():

    def GammaValue(self):
        return((1005*((1013.25)-(0.11986*self._elevation)+(0.000005356*pow(self._elevation,2))))/(622*((2501)-(2.37*self._air_temperature_data))))

class EvapoTranspiration( SoilHeatFlux, NetRadiation, PsychometricConstant, SlopeOfSaturatedVaporPressure):

    def __init__(self, **kwargs):
        if 'Priestley_Taylor_Coefficient' in kwargs: self._priestley_taylor_coefficient = kwargs.get('Priestley_Taylor_Coefficient')
        super().__init__(**kwargs)

    def EvapoTranspirationValue(self):
        net_radiation_value = NetRadiation.NetRadiationValue(self)
        soil_heat_flux_value = SoilHeatFlux.SoilHeatFluxValue(self)
        delta_value = SlopeOfSaturatedVaporPressure.delta(self)
        gamma_value = PsychometricConstant.GammaValue(self)        
        return((self._priestley_taylor_coefficient*(net_radiation_value - soil_heat_flux_value) * delta_value) / (delta_value+gamma_value))

def main():     
    net_radiation_components = EvapoTranspiration(albedo = 10, elevation = 10,air_temperature_data = 27, out_going_long_wave_radiation = 100, insolation = 600, NDVI = 0.7, land_surface_temperature = 45, Priestley_Taylor_Coefficient = 1.26)
    print(net_radiation_components.EvapoTranspirationValue())

if __name__ == '__main__': main()



