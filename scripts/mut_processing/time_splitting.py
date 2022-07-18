from typing import List
import pandas as pd

PATH_TO_MUT = "./data/mutation_dists.filtered.csv"
PATH_TO_METADATA = "./data/gisaid/metadata_clean.csv"
PATH_TO_NEW_MUT = "./data/mutation_dists.filtered.splits.csv"

splits = [
    ["2019-12", "2020-01", "2020-02", "2020-03", "2020-04"],
    ["2020-05", "2020-06", "2020-07"],
    ["2020-08", "2020-09", "2020-10"],
    ["2020-11", "2020-12"],
    ["2021-01", "2021-02"],
    ["2021-03", "2021-04"],
    ["2021-05", "2021-06"],
    ["2021-07", "2021-08"],
    ["2021-09", "2021-10"],
]


def prepare_mut(mut, annot):
    # only terminal mutations
    mut = mut[~mut.child_node.str.startswith("#")]

    # remove ambiguous records from annot table
    _annot_vc = annot.strain.value_counts()
    amb_strains = set(_annot_vc[(_annot_vc > 1)].index)
    annot = annot[~annot.strain.isin(amb_strains)]

    # add columns date and month to mut
    mut2 = pd.merge(
        mut,
        annot[["strain", "date"]],
        left_on="child_node",
        right_on="strain"
    )
    mut2.drop("strain", axis=1, inplace=True)
    mut2 = mut2[mut2.date.str.count("-") > 1]
    mut2["month"] = mut2.date.str.extract("(20..-..)")
    return mut2


def add_splitter_column(mut, splits: List[List[str]]):
    mut = mut.copy()
    mut["date_split"] = -1
    for i, spl in enumerate(splits):
        mut["date_split"][mut.month.isin(set(spl))] = i
    return mut


def main():
    mut = pd.read_csv(PATH_TO_MUT)
    annot = pd.read_csv(PATH_TO_METADATA)

    mut = prepare_mut(mut, annot)
    mut = add_splitter_column(mut, splits)

    mut.to_csv(PATH_TO_NEW_MUT, index=None)


if __name__ == "__main__":
    main()

# month      number_of_mutations
# 2019-12       18
# 2020-01      215
# 2020-02      283
# 2020-03     7241
# 2020-04     7294

# 2020-05     4470
# 2020-06     5076
# 2020-07     6501

# 2020-08     6653
# 2020-09     8345
# 2020-10    12203

# 2020-11    16789
# 2020-12    23739

# 2021-01    41988
# 2021-02    43626

# 2021-03    55907
# 2021-04    59357

# 2021-05    42112
# 2021-06    29576

# 2021-07    54306
# 2021-08    70428

# 2021-09    43771
# 2021-10     2870
