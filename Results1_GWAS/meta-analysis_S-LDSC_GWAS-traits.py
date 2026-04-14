#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob
import os
import re
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def remove_suffix(annotation):
    """Remove L2_* suffix from annotation names"""
    return re.sub(r'L2_\d+$', '', annotation)

def random_effects_meta_analysis(enrichments, std_errors):
    """
    Perform random-effects meta-analysis using DerSimonian-Laird method

    Parameters:
    enrichments: array of enrichment values
    std_errors: array of standard errors

    Returns:
    dict with meta-analysis results
    """
    enrichments = np.array(enrichments)
    std_errors = np.array(std_errors)

    # Calculate weights (inverse variance)
    weights = 1 / (std_errors ** 2)

    # Fixed-effects estimate
    fe_estimate = np.sum(weights * enrichments) / np.sum(weights)
    fe_variance = 1 / np.sum(weights)
    fe_std_error = np.sqrt(fe_variance)

    # Calculate Q statistic for heterogeneity
    Q = np.sum(weights * (enrichments - fe_estimate) ** 2)
    df = len(enrichments) - 1

    # Calculate tau-squared (between-study variance)
    if df > 0:
        tau_squared = max(0, (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights)))
    else:
        tau_squared = 0

    # Random-effects weights
    re_weights = 1 / (std_errors ** 2 + tau_squared)

    # Random-effects estimate
    re_estimate = np.sum(re_weights * enrichments) / np.sum(re_weights)
    re_variance = 1 / np.sum(re_weights)
    re_std_error = np.sqrt(re_variance)

    # Calculate z-score and p-value
    z_score = re_estimate / re_std_error
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    # Calculate I-squared (heterogeneity measure)
    if Q > 0:
        i_squared = max(0, (Q - df) / Q) * 100
    else:
        i_squared = 0

    return {
        'n_studies': len(enrichments),
        'enrichment': re_estimate,
        'std_error': re_std_error,
        'z_score': z_score,
        'p_value': p_value,
        'tau_squared': tau_squared,
        'i_squared': i_squared,
        'q_statistic': Q,
        'q_p_value': 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0
    }

def main():
    print("Starting random-effects meta-analysis...")

    # Step 1: Read significant annotations file
    print("Reading significant annotations file...")
    try:
        annots_df = pd.read_csv('../all_annots_assessment.tsv', sep='\t')
        print(f"Loaded {len(annots_df)} annotations")
    except FileNotFoundError:
        print("Error: Could not find ../all_annots_assessment.tsv")
        return
    except Exception as e:
        print(f"Error reading annotations file: {e}")
        return

    # Step 2: Clean annotation names (remove L2_* suffix)
    clean_annotations = [remove_suffix(ann) for ann in annots_df['Annotation']]
    print(f"Cleaned annotation names: {clean_annotations[:5]}..." if len(clean_annotations) > 5 else f"Cleaned annotation names: {clean_annotations}")

    # Step 3: Find all results files
    results_pattern = 'results_by_trait/*.results'
    results_files = glob.glob(results_pattern, recursive=True)

    if not results_files:
        print(f"No results files found matching pattern: {results_pattern}")
        return

    print(f"Found {len(results_files)} results files")

    # Step 4: Process each results file and collect data
    meta_data = {}  # Dictionary to store data for each annotation

    for results_file in results_files:
        print(f"Processing {results_file}...")

        try:
            results_df = pd.read_csv(results_file, sep='\t')

            # Clean category names (remove L2_* suffix)
            results_df['Clean_Category'] = results_df['Category'].apply(remove_suffix)

            # Filter for categories that match our annotations
            matching_results = results_df[results_df['Clean_Category'].isin(clean_annotations)]

            print(f"  Found {len(matching_results)} matching categories")

            # Store data for meta-analysis
            for _, row in matching_results.iterrows():
                category = row['Clean_Category']
                enrichment = row['Enrichment']
                std_error = row['Enrichment_std_error']

                # Skip if enrichment or std_error is NaN
                if pd.isna(enrichment) or pd.isna(std_error):
                    continue

                if category not in meta_data:
                    meta_data[category] = {'enrichments': [], 'std_errors': [], 'files': []}

                meta_data[category]['enrichments'].append(enrichment)
                meta_data[category]['std_errors'].append(std_error)
                meta_data[category]['files'].append(os.path.basename(os.path.dirname(results_file)))

        except Exception as e:
            print(f"  Error processing {results_file}: {e}")
            continue

    # Step 5: Perform meta-analysis for each annotation
    print("\nPerforming random-effects meta-analysis...")

    meta_results = []

    for annotation in clean_annotations:
        if annotation in meta_data:
            data = meta_data[annotation]

            if len(data['enrichments']) >= 2:  # Need at least 2 studies for meta-analysis
                try:
                    result = random_effects_meta_analysis(data['enrichments'], data['std_errors'])
                    result['annotation'] = annotation
                    result['studies'] = '; '.join(data['files'])
                    meta_results.append(result)

                    print(f"  {annotation}: {result['n_studies']} studies, enrichment = {result['enrichment']:.4f} ± {result['std_error']:.4f}")

                except Exception as e:
                    print(f"  Error in meta-analysis for {annotation}: {e}")
            else:
                print(f"  {annotation}: insufficient data ({len(data['enrichments'])} studies)")
        else:
            print(f"  {annotation}: no matching data found")

    # Step 6: Save results
    if meta_results:
        results_df = pd.DataFrame(meta_results)

        # Reorder columns
        column_order = ['annotation', 'n_studies', 'enrichment', 'std_error', 'z_score', 'p_value',
                       'tau_squared', 'i_squared', 'q_statistic', 'q_p_value', 'studies']

        results_df = results_df[column_order]

        # Sort by p-value
        results_df = results_df.sort_values('p_value')

        # Save to file
        output_file = 'random-effect_meta-analysis.tsv'
        results_df.to_csv(output_file, sep='\t', index=False)

        print(f"\nResults saved to {output_file}")
        print(f"Meta-analysis completed for {len(meta_results)} annotations")

        # Display summary
        print("\nSummary of results:")
        print(results_df[['annotation', 'n_studies', 'enrichment', 'std_error', 'p_value']].head(10))

    else:
        print("No meta-analysis results to save")

if __name__ == "__main__":
    main()
