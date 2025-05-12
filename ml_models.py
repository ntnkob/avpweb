from pycaret.classification import load_model, predict_model
from protlearn.features import aac, ngram, paac, apaac, ctdc, ctdt, ctdd
from helper import reformat_result
import pandas as pd
import numpy as np

def generate_features(idList, seqList):
    def get_features(seq):
        feature = []
        feature_names = []

        feature.extend(aac(seq)[0])
        feature_names.extend(aac(seq)[1])

        feature.extend(ngram(seq)[0])
        feature_names.extend(ngram(seq)[1])

        feature.extend(paac(seq, lambda_ = 5, w = 0.05)[0])
        feature_names.extend(paac(seq, lambda_ = 5, w = 0.05)[1])

        feature.extend(apaac(seq, lambda_ = 5, w = 0.05)[0])
        feature_names.extend(apaac(seq, lambda_ = 5, w = 0.05)[1])

        feature.extend(ctdc(seq)[0])
        feature_names.extend(ctdc(seq)[1])

        feature.extend(ctdt(seq)[0])
        feature_names.extend(ctdt(seq)[1])

        feature.extend(ctdd(seq)[0])
        feature_names.extend(ctdd(seq)[1])
            
        return feature, feature_names
    
    def fix_columns(inp_df):
        passed_columns = []
        for colname in inp_df.columns:
            if colname not in passed_columns:
                passed_columns.append(colname)
            elif f'{colname}.1' not in passed_columns:
                passed_columns.append(f'{colname}.1')
            elif f'{colname}.2' not in passed_columns:
                passed_columns.append(f'{colname}.2')
            else:
                print(f'Too much repetition: {colname}')
        return passed_columns

    all_features = []
    for seq in seqList:
        feature, feature_names = get_features(seq)
        all_features.append(np.concatenate(feature))
    feature_df = pd.DataFrame(np.vstack(all_features), columns = feature_names)
    feature_df.index = idList
    feature_df.columns = fix_columns(feature_df)
    return feature_df

def get_prediction(all_fasta: dict, avp_threshold: float = 0.5, anticov_threshold: float = 0.5, rename_dict: dict = {}) -> None: 
    
    feature_df = generate_features(all_fasta['id'], all_fasta['seq'])

    # Step 1: AVP vs Non-AVP
    AVP_model = load_model('assets/models/avp_model')
    AVP_result = predict_model(AVP_model, data = feature_df, probability_threshold=avp_threshold)
    importance_sorted_arg = AVP_model._final_estimator.feature_importances_.argsort()
    AVP_result = AVP_result.iloc[:, list(importance_sorted_arg) + [748, 749]]
    AVP_reformat = reformat_result(AVP_result, rename_dict = {'SeqID':'Sequence ID', 'prediction_label': 'Class', 'prediction_score':'Probability'})
    AVP_reformat.to_csv('output/avp_result.csv', index=False)

    # Step 2: Filter only AVPs
    AVP_df = AVP_result[AVP_result['prediction_label']==1].drop(columns = ['prediction_label','prediction_score'])
    AntiCoV_model = load_model('assets/models/anticov_model')
    AntiCoV_result = predict_model(AntiCoV_model, data = AVP_df, probability_threshold=anticov_threshold)
    importance_sorted_arg = AntiCoV_model._final_estimator.feature_importances_.argsort()
    AntiCoV_result = AntiCoV_result.iloc[:, list(importance_sorted_arg) + [748, 749]]
    AntiCoV_reformat = reformat_result(AntiCoV_result, rename_dict = {'SeqID':'Sequence ID', 'prediction_label': 'Class', 'prediction_score':'Probability'})
    AntiCoV_reformat.to_csv('output/anticov_result.csv', index=False)