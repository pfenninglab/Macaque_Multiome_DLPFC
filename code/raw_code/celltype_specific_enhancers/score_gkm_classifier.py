import os, sys, gc, math, argparse
from sklearn import metrics
import pandas as pd
import numpy as np

np.seterr(divide = 'ignore')

def read_predictions(pos, neg):
    pos_df = pd.read_table(pos, names= ['ids', 'y_pred_score']); pos_df['y'] = 1
    neg_df = pd.read_table(neg, names= ['ids', 'y_pred_score']); neg_df['y'] = 0
    df = pd.concat([pos_df, neg_df]).fillna(0)
    df['y_pred_score'] = np.nan_to_num(df['y_pred_score'])
    df['y_pred_class'] = df['y_pred_score'] > 0
    return df

def evaluate_sequences(model_name, pos, neg, args):
    # Creating an empty Dataframe with column names only
    df = pd.DataFrame(vars(args), index=[0])
    tmp = read_predictions(pos, neg)
    # compute prediction statistics 
    tn, fp, fn, tp = metrics.confusion_matrix(tmp['y']==1, tmp['y_pred_class']).ravel()
    accuracy = metrics.balanced_accuracy_score(tmp['y']==1, tmp['y_pred_class'])
    f1_score = metrics.f1_score(tmp['y']==1, tmp['y_pred_class'], average = 'weighted')
    fhalf_score = metrics.fbeta_score(tmp['y']==1, tmp['y_pred_class'], beta = 0.5, average = 'weighted')
    roc_auc = metrics.roc_auc_score(tmp['y']==1, tmp['y_pred_score'])
    precision, recall, thresholds = metrics.precision_recall_curve(tmp['y']==1, tmp['y_pred_score'])
    prc_auc = metrics.auc(recall, precision) # x, y
    print(f'Accuracy: {accuracy}. ')
    print(f'f1_score: {f1_score}.')
    print(f'roc_auc: {roc_auc}.')
    print(f'prc_auc: {prc_auc}.')
    # add this row to dataframe
    df = pd.concat([df.reset_index(drop=True),
        pd.DataFrame({'model': model_name, 
            'tn' : tn, 'fp' : fp, 
            'fn': fn, 'tp' : tp,
            'accuracy': accuracy, 
            'auROC': roc_auc, 
            'auPRC': prc_auc,
            'f1_score': f1_score, 
            'fhalf_score': fhalf_score}, index=[0])], axis = 1)    
    return df

def main(args):
    """Main function
    Args:
    args (argparse):
    """
    # call main functions
    print('In evaluation mode.')
    if not os.path.exists(args.model_name):
        print('No model found with specified training parameters. Please train model first.')
        return
    predict_out = os.path.basename(args.model_name).replace('.model.txt', '')
    model_test_performance = f'{args.out_dir}/predictions/{args.prefix}/{predict_out}_gkmpredict_{args.label}_eval.feather'
    print(f'Model performance to be written to {model_test_performance}')
    if os.path.exists(model_test_performance) and not args.force:
        print(f'The performance file exists w/o permission to overwrite. Use --force to overwrite.')
        return
    df = evaluate_sequences(args.model_name, args.pred_pos, args.pred_neg, args)
    # save model performances to feather object
    if not os.path.exists(f'{args.out_dir}/predictions/{args.prefix}'):
        os.makedirs(f'{args.out_dir}/predictions/{args.prefix}')
    df.to_feather(model_test_performance)
    return


if __name__ == '__main__':  
    #### set cnn parameters:
    parser = argparse.ArgumentParser(description='Parse GKM model evaluation parameters.')
    parser.add_argument("--prefix", type=str, help="sub folder to put predictions inside outdir")
    parser.add_argument("--model_name", type=str, help="complete model name")
    parser.add_argument("--label", type=str, help="discriminative label of the prediciton files, e.g. valid, test, other files, etc.")
    parser.add_argument("--pred_pos", type=str, help="prediction file of positives.")
    parser.add_argument("--pred_neg", type=str, help="prediction file of negatives.")
    parser.add_argument("--out_dir", type=str, default = '.', help="path to ouputput directory, default is pwd")
    parser.add_argument("--force", help="Whether to overwrite previously scored files.", action='store_true')

    ### parse arguments
    args = parser.parse_args(['--prefix=L2.CUX2.MEIS2', '--label=valid',
        '--out_dir=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/celltype_specific_enhancers', 
        '--model_name=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/celltype_specific_enhancers/models_svm/L2.CUX2.MEIS2/L2.CUX2.MEIS2vsEXC_fold1_t4_l7_k6_d1_c5_w1.5.model.txt',
        '--pred_pos=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/celltype_specific_enhancers/predictions/L2.CUX2.MEIS2/L2.CUX2.MEIS2vsEXC_fold1_t4_l7_k6_d1_c5_w1.5_gkmpredict_valid_positive.txt', 
        '--pred_neg=/projects/pfenninggroup/singleCell/Macaque_Multiome_DLPFC/data/tidy_data/celltype_specific_enhancers/predictions/L2.CUX2.MEIS2/L2.CUX2.MEIS2vsEXC_fold1_t4_l7_k6_d1_c5_w1.5_gkmpredict_valid_negative.txt'])

    args = parser.parse_args()
    main(args)