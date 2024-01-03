import swifter


def split_column_into_columns(table, column_name, delimeter=';'):
    # Build mutliple copies for inter protein results
    # https://blog.csdn.net/Elimeny/article/details/93030861
    target_column = table[column_name].str.split(delimeter, expand=True)
    target_column = target_column.stack()
    target_column = target_column.reset_index(
        level=1, drop=True, name=column_name)
    target_column = target_column.rename(column_name)
    table = table.drop([column_name], axis=1).join(target_column)

    # appended_data = pd.DataFrame()
    # for _, line in tqdm.tqdm(PSM_res[PSM_res[column_name].swifter.apply(lambda x: len(re.split(";", x))) > 1].iterrows()):
    #     for i, protein in enumerate(re.split(";", line[column_name])):
    #         if i == 0:
    #             continue
    #         new_line = line.copy()
    #         new_line.at[column_name] = protein
    #         appended_data = appended_data.append(new_line, ignore_index=True)
    # PSM_res[column_name] = PSM_res[column_name].swifter.apply(
    #     lambda x: re.split(";", x)[0] if len(re.split(";", x)) > 1 else x)

    # # Merge copies together
    # PSM_res = PSM_res.append(appended_data, ignore_index=True)
    return table
    # TODO: split data and split columns


def group_by_col(df, by, keep_min_cols=["Evalue", "Score"],
                 keep_max_cols=[], sum_cols=[],
                 concat_str_cols=[]):

    def merge_rule(x, keep_min_cols, keep_max_cols=[], sum_cols=[], concat_str_cols=[]):
        if len(x) == 1:
            return x.iloc[0]

        if not keep_min_cols == []:
            df = x.sort_values(by=keep_min_cols, ascending=True)
        elif not keep_max_cols == []:
            df = x.sort_values(by=keep_max_cols, ascending=False)
        else:
            df = x
        for col in keep_max_cols:
            df[col].iloc[0] = max(df[col])
        for col in sum_cols:
            df[col].iloc[0] = sum(df[col])
        for col in concat_str_cols:
            df[col].iloc[0] = ";".join(df[col])
        return df.iloc[0]

    return df.groupby(by, dropna=False).apply(
        merge_rule, keep_min_cols, keep_max_cols, sum_cols, concat_str_cols)
    return df.groupby(by, dropna=False).agg({
        concat_str_cols: lambda x: ';'.join(x),
        keep_min_cols: 'min',
        keep_max_cols: 'max',
        sum_cols: 'sum',
        })


def group_by_col_swifter(df, by, keep_min_cols=["Evalue", "Score"],
                         keep_max_cols=[], sum_cols=[],
                         concat_str_cols=[]):

    def merge_rule(x, keep_min_cols, keep_max_cols=[], sum_cols=[], concat_str_cols=[]):
        if len(x) == 1:
            return x.iloc[0]

        if not keep_min_cols == []:
            df = x.sort_values(by=keep_min_cols, ascending=True)
        elif not keep_max_cols == []:
            df = x.sort_values(by=keep_max_cols, ascending=False)
        else:
            df = x
        for col in keep_max_cols:
            df[col].iloc[0] = max(df[col])
        for col in sum_cols:
            df[col].iloc[0] = sum(df[col])
        for col in concat_str_cols:
            df[col].iloc[0] = ";".join(df[col])
        return df.iloc[0]

    return df.swifter.groupby(by, dropna=False).apply(
        merge_rule, keep_min_cols, keep_max_cols, sum_cols, concat_str_cols)
