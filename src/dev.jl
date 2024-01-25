using IdealGases

gas1d = IdealGases.Gas1D()
gas = Gas()

gas2 = IdealGases.fuel_combustion(gas1d, "H2", 298.15, 0.01260625127079659)

gas3 = IdealGases.gas_mixing(gas, gas2, 2.0)