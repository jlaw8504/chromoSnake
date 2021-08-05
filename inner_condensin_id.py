from pickle import load
import pandas as pd

with open('/home/test/PycharmProjects/chromoSnake/equib_df_sizes.pkl', 'rb') as f:
    df_sizes = load(f)

# append condensin type column to data frame
df_sizes.loc[:, 'condensin_type'] = pd.Series([None] * df_sizes.shape[0], index=df_sizes.index)

chromosome_list = [f'chr{i}_' for i in range(1, 17)]
for chromosome in chromosome_list:
    for loop_size in pd.unique(df_sizes.loop_size):
        df_temp = df_sizes.loc[
            (df_sizes.loop_size == loop_size) &
            (
                (
                    (df_sizes.mass_label_1.str.match('super')) &
                    (df_sizes.mass_label_2.str.contains(chromosome))
                ) |
                (
                    (df_sizes.mass_label_1.str.contains(chromosome)) &
                    (df_sizes.mass_label_2.str.match('super'))
                )
            )
        ]
        mass_index_list = list(df_temp.spring_mass_idx_1) + list(df_temp.spring_mass_idx_2)
        mass_index_list.sort()
        strand_1_range = range(mass_index_list[0], mass_index_list[3]+1)
        strand_2_range = range(mass_index_list[4], mass_index_list[7]+1)


