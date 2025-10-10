import CoolProp.CoolProp as CP

P0 = 692*133.322
T0 = 24+273.15
work_fluid = 'air'
air_rho = CP.PropsSI('D', 'T', T0, 'P', P0, work_fluid)
data={ # Frequency: [[Q, H]]
    '30Hz': [
        [49.717139035644166, 74.5],   
        [78.12830001633337, 71],      
        [118.75128159554802, 55.5],   
        [126.59662786498599, 52.5],   
        [100.47926032217654, 60]      
    ],
    '60Hz': [
        [104.10988439084795, 297.5],  
        [150.58516712918168, 278],    
        [182.43990133394215, 259],    
        [216.68391878025795, 236],    
        [245.85917191003233, 214]     
    ]
}

for e in data:
    for num in range(0, len(data[e])):
        data[e][num][1]*=1/1000 # mmH20 to mH20

for value in data.values():
    print(value)
