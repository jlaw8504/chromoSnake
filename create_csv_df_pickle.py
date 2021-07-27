from tensionSim import TensionSim as ts
import pandas as pd
from pickle import dump
import datetime as dt

start = dt.datetime.now()
my_tension = ts('/home/test/Documents/chromoShake/varSize/WTspindle_2_25_coh220.out')
my_tension.parse_sim('/home/test/PycharmProjects/chromoSnake/new_test.txt')
names_list = [
    'spring_idx',
    'time',
    'spring_mass_idx_2',
    'spring_mass_idx_1',
    'tension_N',
    'mass_label_1',
    'mass_label_2',
    'is_in_loop_1',
    'is_in_loop_2'
    ]
df = pd.read_csv('/home/test/PycharmProjects/chromoSnake/new_test.txt', names=names_list)
equilibrium_df = df.loc[df['time'] > 0.02]
del df
calc_df = equilibrium_df[['spring_idx', 'tension_N']]
mean_calc_df = calc_df.groupby(['spring_idx']).mean()
label_df = equilibrium_df.loc[equilibrium_df['time'] == equilibrium_df['time'].iloc[0]]
del equilibrium_df
label_df = label_df.set_index(['spring_idx'])
label_df.drop(columns=['time', 'tension_N'], inplace=True)
final_df = label_df.join(mean_calc_df)
with open('/home/test/PycharmProjects/chromoSnake/WTspindle_2_25_coh220_final_df.pkl', 'wb') as f_pkl:
    dump(final_df, f_pkl)

end = dt.datetime.now()
print((end-start))
