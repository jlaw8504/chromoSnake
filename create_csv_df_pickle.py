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


# add condensin labels
def id_condensin(df):
    if abs(df['spring_mass_idx_1'] - df['spring_mass_idx_2']) != 1 and 'cohesin' not in df['mass_label_1']:
        return True
    else:
        return False


final_df['is_condensin_spring'] = final_df.apply(id_condensin, axis=1)
# convert the tension from N to pN
final_df['tension_pN'] = final_df.apply(lambda x: x['tension_N']*10**12, axis=1)
final_df.drop(columns=['tension_N'], inplace=True)

with open('/home/test/PycharmProjects/chromoSnake/WTspindle_2_25_coh220_final_df.pkl', 'wb') as f_pkl:
    dump(final_df, f_pkl)

end = dt.datetime.now()
print((end - start))
