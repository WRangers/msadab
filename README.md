# Codes of *Modelling and Simulation Analysis of Drones Allocation for Bushfires*

## File Tree

```
msadba
│  README.md
│  
├─code
│      clean_data.m				# Data preprocessing script
|      cruise_line.m             # Script file to design the cruise lines
│      data_proces_2019.m        # Process script file of data of year 2019
│      demo.m                    # Demo script to show simulation
│      mycmap.mat                # The custom color map data
│      
└─data
        blo_fire_data.mat        # Processed data
        modis_2019_Australia.csv # Bushfire data collected by satellites
        victoria.csv             # The border data of Victoria, Australia
```

## Testing the algorithm

1. Load the data `blo_fire_data.mat `.
2. Run the script `demo.m`.

## Adapting to other data

1. Change the data `victoria.csv` to the desire region's.
2. Run the script `clean_data.m` to clean the data.
3. Run the script `data_proces_2019.m` to get the data needed. Attention: you may come across some error due the data format, please keep the same format with that of data `modis_2019_Australia.csv`'s.
4. Run the script `demo.m`.