import pandas as pd

def reformat_result(input_df: pd.DataFrame, rename_dict: dict = {}) -> pd.DataFrame:
    '''
    Rename and reorder the columns in `input_df` according to `rename_dict`
    '''
    final_result = input_df.reset_index(names='SeqID')
    column_list = final_result.columns.tolist()
    for location, original_name in enumerate(rename_dict):
        column_list.insert(location, column_list.pop(column_list.index(original_name)))
    final_result = final_result.loc[:, column_list]
    final_result = final_result.rename(columns = rename_dict)
    final_result['Class'] = final_result['Class'].replace({0: "Non-AVP", 1: "AVP"})
    return final_result