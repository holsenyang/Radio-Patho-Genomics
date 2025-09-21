from sklearn.metrics import roc_curve, roc_auc_score, accuracy_score, recall_score, precision_score, f1_score
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import chi2
from sklearn.metrics import confusion_matrix
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

def bootstrap_auc(y_true, y_pred, n_bootstrap=2000, alpha=0.05, random_seed=0):
    np.random.seed(random_seed)
    auc_values = []
    n = len(y_true)
    for _ in range(n_bootstrap):
        indices = np.random.choice(n, size=n, replace=True)
        y_true_bootstrap = y_true[indices]
        y_pred_bootstrap = y_pred[indices]
        auc = roc_auc_score(y_true_bootstrap, y_pred_bootstrap)
        auc_values.append(auc)
    auc_values = np.array(auc_values)
    lower = np.percentile(auc_values, 100 * alpha / 2)
    upper = np.percentile(auc_values, 100 * (1 - alpha / 2))
    return lower, upper

def delong_test(y_true, pred1, pred2):
    auc1 = roc_auc_score(y_true, pred1)
    auc2 = roc_auc_score(y_true, pred2)
    
    fpr1, tpr1, _ = roc_curve(y_true, pred1)
    fpr2, tpr2, _ = roc_curve(y_true, pred2)
    
    def calculate_covariance(y_true, pred1, pred2):
        n = len(y_true)
        sigma = np.zeros((2, 2))
        for i in range(n):
            for j in range(n):
                if y_true[i] == 1 and y_true[j] == 1:
                    sigma[0, 0] += (pred1[i] - pred1[j]) * (pred1[i] - pred1[j])
                    sigma[0, 1] += (pred1[i] - pred1[j]) * (pred2[i] - pred2[j])
                    sigma[1, 1] += (pred2[i] - pred2[j]) * (pred2[i] - pred2[j])
        sigma /= n * (n - 1)
        return sigma
    
    sigma = calculate_covariance(y_true, pred1, pred2)
    z = (auc1 - auc2) / np.sqrt(sigma[0, 0] + sigma[1, 1] - 2 * sigma[0, 1])
    p_value = 2 * (1 - norm.cdf(abs(z)))
    return p_value

def find_matching_threshold(y_true, y_pred, target_specificity=None, target_sensitivity=None):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    
    if target_specificity is not None:
        specificity = 1 - fpr
        abs_diff = np.abs(specificity - target_specificity)
        optimal_idx = np.argmin(abs_diff)
        optimal_threshold = thresholds[optimal_idx]
        matched_sensitivity = tpr[optimal_idx]
        return optimal_threshold, matched_sensitivity
    
    if target_sensitivity is not None:
        abs_diff = np.abs(tpr - target_sensitivity)
        optimal_idx = np.argmin(abs_diff)
        optimal_threshold = thresholds[optimal_idx]
        matched_specificity = 1 - fpr[optimal_idx]
        return optimal_threshold, matched_specificity

def bootstrap_ci(y_true, y_pred, metric_func, n_bootstrap=2000, alpha=0.05):
    np.random.seed(0)
    metric_values = []
    n = len(y_true)
    for _ in range(n_bootstrap):
        indices = np.random.choice(n, size=n, replace=True)
        y_true_bootstrap = y_true[indices]
        y_pred_bootstrap = y_pred[indices]
        metric_value = metric_func(y_true_bootstrap, y_pred_bootstrap)
        metric_values.append(metric_value)
    metric_values = np.array(metric_values)
    lower = np.percentile(metric_values, 100 * alpha / 2)
    upper = np.percentile(metric_values, 100 * (1 - alpha / 2))
    return lower, upper

def mcnemar_test(model_pred, doctor_pred, gt):
    gt_pos_mask = (gt == 1)
    model_pred_pos = model_pred[gt_pos_mask]
    doctor_pred_pos = doctor_pred[gt_pos_mask]
    sensitive_cm = confusion_matrix(doctor_pred_pos, model_pred_pos)

    gt_neg_mask = (gt == 0)
    model_pred_neg = model_pred[gt_neg_mask]
    doctor_pred_neg = doctor_pred[gt_neg_mask]
    specific_cm = confusion_matrix(doctor_pred_neg, model_pred_neg)

    def mcnemar_p(cm):
        b = cm[0, 1]
        c = cm[1, 0]
        if (b + c) == 0:
            return 1.0
        stat = (b - c)**2 / (b + c)
        return 1 - chi2.cdf(stat, df=1)

    sensitive_p = mcnemar_p(sensitive_cm)
    specific_p = mcnemar_p(specific_cm)
    return sensitive_p, specific_p

def compare_ai_vs_doctors(my_data, model_names_all, task, subset='AI VS Pathologist'):
    ALL_results = my_data[my_data['Group'].isin(['test'])]
    model_names = model_names_all
    pred_column = [f'{task}-0', f'{task}-1']

    gt = np.array(ALL_results[task])
    preds = [np.array(ALL_results[d]) for d in model_names]

    plt.figure(figsize=(6, 5))

    metric_list = []
    for i, (mname, pred) in enumerate(zip(model_names, preds)):
        if set(np.unique(pred)) == {0, 1}:
            sensitivity = recall_score(gt, pred)
            specificity = recall_score(1 - gt, 1 - pred)
            accuracy = accuracy_score(gt, pred)
            ppv = precision_score(gt, pred)
            npv = precision_score(1 - gt, 1 - pred)
            precision = precision_score(gt, pred)
            recall = recall_score(gt, pred)
            f1 = f1_score(gt, pred)
            threshold = 0.5
            
            if i == 1:
                plt.scatter(1 - specificity, sensitivity, marker='^', color='green', label=f'{mname} (Sensitivity: {sensitivity:.3f}, Specificity: {specificity:.3f})')
            elif i == 2:
                plt.scatter(1 - specificity, sensitivity, marker='s', color='purple', label=f'{mname} (Sensitivity: {sensitivity:.3f}, Specificity: {specificity:.3f})')
            
            auc_value = roc_auc_score(gt, pred)
            ci_lower, ci_upper = bootstrap_auc(gt, pred)
            ci = f"{ci_lower:.3f} - {ci_upper:.3f}"
            
            metric_list.append([mname, accuracy, auc_value, ci, sensitivity, specificity, ppv, npv, precision, recall, f1, threshold, subset])
            
        else:
            fpr, tpr, thresholds = roc_curve(gt, pred)
            auc_value = roc_auc_score(gt, pred)
            ci_lower, ci_upper = bootstrap_auc(gt, pred)
            ci = f"{ci_lower:.3f} - {ci_upper:.3f}"
            
            plt.plot(fpr, tpr, label=f'{mname} AUC: {auc_value:.3f} ({ci})', linewidth=2)
            
            youden_index = tpr - fpr
            optimal_idx = np.argmax(youden_index)
            optimal_threshold = thresholds[optimal_idx]
            pred_labels = np.where(pred >= optimal_threshold, 1, 0)
            
            accuracy = accuracy_score(gt, pred_labels)
            sensitivity = tpr[optimal_idx]
            specificity = 1 - fpr[optimal_idx]
            ppv = precision_score(gt, pred_labels)
            npv = precision_score(1 - gt, 1 - pred_labels)
            precision = precision_score(gt, pred_labels)
            recall = recall_score(gt, pred_labels)
            f1 = f1_score(gt, pred_labels)
            
            metric_list.append([mname, accuracy, auc_value, ci, sensitivity, specificity, ppv, npv, precision, recall, f1, optimal_threshold, subset])

    has_doctors = any(set(np.unique(pred)) == {0, 1} for pred in preds[1:])

    if has_doctors:
        human_preds = [pred for pred in preds[1:] if set(np.unique(pred)) == {0, 1}]
        human_aucs = []
        for pred in human_preds:
            auc = roc_auc_score(gt, pred)
            human_aucs.append(auc)
        average_auc = np.mean(human_aucs)

        def bootstrap_mean_auc(aucs, n_bootstrap=2000, alpha=0.05):
            mean_aucs = []
            for _ in range(n_bootstrap):
                indices = np.random.choice(len(aucs), size=len(aucs), replace=True)
                bootstrap_aucs = np.array(aucs)[indices]
                mean_auc = np.mean(bootstrap_aucs)
                mean_aucs.append(mean_auc)
            mean_aucs = np.array(mean_aucs)
            lower = np.percentile(mean_aucs, 100 * alpha / 2)
            upper = np.percentile(mean_aucs, 100 * (1 - alpha / 2))
            return lower, upper

        human_ci_lower, human_ci_upper = bootstrap_mean_auc(human_aucs)
        human_ci = f"{human_ci_lower:.3f} - {human_ci_upper:.3f}"

        human_sensitivities = []
        human_specificities = []
        for pred in human_preds:
            sensitivity = recall_score(gt, pred)
            specificity = recall_score(1 - gt, 1 - pred)
            human_sensitivities.append(sensitivity)
            human_specificities.append(specificity)
        average_sensitivity = np.mean(human_sensitivities)
        average_specificity = np.mean(human_specificities)

        def bootstrap_mean_sensitivity_specificity(sensitivities, specificities, n_bootstrap=2000, alpha=0.05):
            mean_sensitivities = []
            mean_specificities = []
            for _ in range(n_bootstrap):
                indices = np.random.choice(len(sensitivities), size=len(sensitivities), replace=True)
                bootstrap_sensitivities = np.array(sensitivities)[indices]
                bootstrap_specificities = np.array(specificities)[indices]
                mean_sensitivity = np.mean(bootstrap_sensitivities)
                mean_specificity = np.mean(bootstrap_specificities)
                mean_sensitivities.append(mean_sensitivity)
                mean_specificities.append(mean_specificity)
            mean_sensitivities = np.array(mean_sensitivities)
            mean_specificities = np.array(mean_specificities)
            lower_sensitivity = np.percentile(mean_sensitivities, 100 * alpha / 2)
            upper_sensitivity = np.percentile(mean_sensitivities, 100 * (1 - alpha / 2))
            lower_specificity = np.percentile(mean_specificities, 100 * alpha / 2)
            upper_specificity = np.percentile(mean_specificities, 100 * (1 - alpha / 2))
            return lower_sensitivity, upper_sensitivity, lower_specificity, upper_specificity

        ci_sensitivity_lower, ci_sensitivity_upper, ci_specificity_lower, ci_specificity_upper = bootstrap_mean_sensitivity_specificity(human_sensitivities, human_specificities)

        plt.errorbar(1 - average_specificity, average_sensitivity, 
                     xerr=[[average_specificity - ci_specificity_lower], [ci_specificity_upper - average_specificity]], 
                     yerr=[[average_sensitivity - ci_sensitivity_lower], [ci_sensitivity_upper - average_sensitivity]], 
                     fmt='o', color='orange', ecolor='orange', capsize=2.5, 
                     label=f'Average Human AUC: {average_auc:.3f} ({human_ci})')

        model_pred = preds[0]
        doctors_pred = [pred for pred in preds[1:] if set(np.unique(pred)) == {0, 1}]

        fpr, tpr, thresholds = roc_curve(gt, model_pred)
        youden_index = tpr - fpr
        optimal_idx = np.argmax(youden_index)
        optimal_threshold = thresholds[optimal_idx]
        model_pred_labels = np.where(model_pred >= optimal_threshold, 1, 0)
        doctors_pred_labels = doctors_pred

        comparison_results = []
        for i, doctor_pred_labels in enumerate(doctors_pred_labels):
            doctor_name = model_names[i + 1]
            
            doctor_sensitivity = recall_score(gt, doctor_pred_labels)
            doctor_specificity = recall_score(1 - gt, 1 - doctor_pred_labels)
            
            model_threshold, model_sensitivity = find_matching_threshold(gt, model_pred, target_specificity=doctor_specificity)
            model_pred_labels_sensitivity = np.where(model_pred >= model_threshold, 1, 0)
            model_sensitivity_actual = recall_score(gt, model_pred_labels_sensitivity)
            ci_lower_sensitivity, ci_upper_sensitivity = bootstrap_ci(gt, model_pred_labels_sensitivity, recall_score)
            
            model_threshold, model_specificity = find_matching_threshold(gt, model_pred, target_sensitivity=doctor_sensitivity)
            model_pred_labels_specificity = np.where(model_pred >= model_threshold, 1, 0)
            model_specificity_actual = recall_score(1 - gt, 1 - model_pred_labels_specificity)
            ci_lower_specificity, ci_upper_specificity = bootstrap_ci(gt, model_pred_labels_specificity, lambda y_true, y_pred: recall_score(1 - y_true, 1 - y_pred))
            
            sensitive_p = mcnemar_test(model_pred_labels_sensitivity, doctor_pred_labels, gt)[0]
            specific_p = mcnemar_test(model_pred_labels_specificity, doctor_pred_labels, gt)[1]
            
            comparison_results.append([
                doctor_name,
                doctor_sensitivity,
                f"{model_sensitivity_actual:.3f} ({ci_lower_sensitivity:.3f}-{ci_upper_sensitivity:.3f})",
                sensitive_p,
                doctor_specificity,
                f"{model_specificity_actual:.3f} ({ci_lower_specificity:.3f}-{ci_upper_specificity:.3f})",
                specific_p
            ])
        
        comparison_df = pd.DataFrame(comparison_results, columns=[
            'Doctor_ID', 'Doctor_Sensitivity', 'Model_Sensitivity_95CI', 'Sensitivity_p_value',
            'Doctor_Specificity', 'Model_Specificity_95CI', 'Specificity_p_value'
        ])
        print(comparison_df)

    else:
        pass

    plt.plot([0, 1], [0, 1], 'k--', linewidth=1.5, alpha=0.7, label='_nolegend_')
    plt.title(f'{subset}', fontsize=20, pad=20, fontweight='bold')
    plt.xlabel('False Positive Rate', fontsize=18)
    plt.ylabel('True Positive Rate', fontsize=18)

    ax = plt.gca()
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)

    handles, labels = ax.get_legend_handles_labels()

    if has_doctors:
        new_handles = [handles[0], handles[-1]] + handles[1:-1]
        new_labels = [labels[0], labels[-1]] + labels[1:-1]
    else:
        new_handles = handles
        new_labels = labels

    ax.legend(new_handles, new_labels, loc='lower right', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7, linewidth=1.5)
    plt.xlim([-0.04, 1.04])
    plt.ylim([-0.04, 1.04])
    plt.show()

    metric_df = pd.DataFrame(metric_list, columns=['Signature', 'Accuracy', 'AUC', '95_CI', 'Sensitivity', 'Specificity', 'PPV', 'NPV', 'Precision', 'Recall', 'F1', 'Threshold', 'Cohort'])
    return metric_df, comparison_df if has_doctors else None

if __name__ == "__main__":
    # Example usage
    # my_data = pd.read_csv('your_data.csv')
    # model_names_all = ['model1', 'doctor1', 'doctor2']
    # task = 'cancer_diagnosis'
    # metric_df, comparison_df = compare_ai_vs_doctors(my_data, model_names_all, task)
    pass