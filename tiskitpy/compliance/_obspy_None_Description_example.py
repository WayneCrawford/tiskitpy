from obspy.core.inventory import Inventory, Network, Station, Channel, Response


# MAKE INVENTORY
seisresp = Response.from_paz([], [], 2000, input_units='m/s', output_units='count')
print(seisresp)
print(seisresp.instrument_sensitivity)
print(seisresp.instrument_sensitivity.input_units_description)
print(seisresp.instrument_sensitivity.input_units_description=='None')
print(seisresp.instrument_sensitivity.input_units_description is None)
print(seisresp.response_stages[0])
print(seisresp.response_stages[0].input_units_description=='None')
print(seisresp.response_stages[0].input_units_description is None)

channels=[Channel('LHZ', '00', 0, 0, 0, 0, response=seisresp, dip = -90)]
stations = [Station('STA', 0, 0, 0, channels=channels)]
networks = [Network('XX', stations=stations)]
inv = Inventory(networks=networks)
inv.write('test.inv.xml', 'STATIONXML')