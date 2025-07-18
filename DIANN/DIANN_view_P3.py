#!/usr/bin/env python3
"""
DIA-NN Proteomics Analysis Pipeline
Includes filtering, statistical analysis, visualization, and GO enrichment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import sys

# Try to import additional packages
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("Note: scikit-learn not installed. PCA plots will be skipped.")

try:
    import plotly.express as px
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Note: plotly not installed. Interactive plots will be skipped.")

class DIANNAnalyzer:
    def __init__(self, quant_file, metadata_file):
        """Initialize with quantification and metadata files"""
        self.quant_file = quant_file
        self.metadata_file = metadata_file
        self.output_dir = Path(quant_file).parent / "analysis_results"
        self.output_dir.mkdir(exist_ok=True)
        
        print("=" * 80)
        print("DIA-NN PROTEOMICS ANALYSIS PIPELINE")
        print("=" * 80)
        
    def load_data(self):
        """Load quantification data and metadata"""
        print("\n1. Loading data...")
        
        # Load protein quantification
        self.proteins = pd.read_csv(self.quant_file, index_col=0)
        print(f"Loaded {len(self.proteins)} proteins")
        
        # Load metadata
        self.metadata = pd.read_csv(self.metadata_file)
        print(f"Loaded metadata for {len(self.metadata)} samples")
        
        # Get sample columns (exclude annotation columns)
        self.sample_cols = [col for col in self.proteins.columns 
                           if col not in ['Protein.Names', 'Genes', 'N_Samples', 
                                         'Mean_Intensity', 'CV_%']]
        
        # Match samples with metadata
        self._match_samples()
        
    def _match_samples(self):
        """Match sample columns with metadata"""
        print("\nMatching samples with metadata...")
        
        # Debug: Show what we're trying to match
        print(f"Sample columns in protein data: {self.sample_cols[:5]}...")
        print(f"Sample names in metadata: {list(self.metadata['Sample_Name'].head())}")
        
        # Create mapping from sample names to metadata
        self.sample_to_group = {}
        self.sample_to_replicate = {}
        
        for _, row in self.metadata.iterrows():
            sample_name = row['Sample_Name']
            
            # Try different matching strategies
            for col in self.sample_cols:
                # Direct match
                if col == sample_name:
                    self.sample_to_group[col] = row['Group']
                    self.sample_to_replicate[col] = row['Replicate']
                    break
                # Partial match - if sample name is contained in column name
                elif sample_name in col:
                    self.sample_to_group[col] = row['Group']
                    self.sample_to_replicate[col] = row['Replicate']
                    break
                # Try matching just the base name without extensions
                elif sample_name.split('.')[0] in col or col.split('.')[0] in sample_name:
                    self.sample_to_group[col] = row['Group']
                    self.sample_to_replicate[col] = row['Replicate']
                    break
        
        # If no matches found, try to match based on 'Run' column pattern
        if len(self.sample_to_group) == 0:
            print("\nTrying to match based on Run column patterns...")
            for _, row in self.metadata.iterrows():
                run_name = row['Run']
                # Extract just the filename part from the Run path
                run_base = Path(run_name).stem
                
                for col in self.sample_cols:
                    if run_base in col or col in run_base:
                        self.sample_to_group[col] = row['Group']
                        self.sample_to_replicate[col] = row['Replicate']
                        print(f"  Matched: {col} -> {row['Group']}")
                        break
        
        print(f"\nMatched {len(self.sample_to_group)} samples to groups")
        
        if len(self.sample_to_group) == 0:
            print("\nWARNING: No samples could be matched!")
            print("Please check that sample names in metadata match column names in protein data")
            print("\nExample column names from protein data:")
            for col in self.sample_cols[:3]:
                print(f"  - {col}")
            print("\nExample entries from metadata:")
            print(self.metadata[['Run', 'Sample_Name', 'Group']].head(3))
            
        # Get unique groups
        self.groups = sorted(set(self.sample_to_group.values()))
        if self.groups:
            print(f"Groups found: {', '.join(self.groups)}")
        
    def filter_proteins(self):
        """Filter out contaminants and proteins not in all conditions"""
        print("\n2. Filtering proteins...")
        initial_count = len(self.proteins)
        
        # Filter out contaminants (genes starting with ##)
        if 'Genes' in self.proteins.columns:
            mask_contaminants = self.proteins['Genes'].fillna('').str.startswith('##')
            self.proteins = self.proteins[~mask_contaminants]
            print(f"Removed {mask_contaminants.sum()} contaminant proteins")
        
        # Filter proteins that aren't in at least one replicate of every condition
        print("\nFiltering for proteins present in all conditions...")
        
        # Create a matrix of presence/absence for each group
        presence_matrix = pd.DataFrame(index=self.proteins.index, columns=self.groups)
        
        for group in self.groups:
            # Get samples for this group
            group_samples = [col for col, g in self.sample_to_group.items() if g == group]
            # Check if protein is present in at least one replicate
            presence_matrix[group] = (~self.proteins[group_samples].isna()).sum(axis=1) >= 1
        
        # Keep only proteins present in all groups
        mask_all_groups = presence_matrix.all(axis=1)
        self.proteins = self.proteins[mask_all_groups]
        
        print(f"Removed {(~mask_all_groups).sum()} proteins not in all conditions")
        print(f"Final protein count: {len(self.proteins)} (removed {initial_count - len(self.proteins)} total)")
        
        # Save filtered data
        output_file = self.output_dir / "filtered_proteins.csv"
        self.proteins.to_csv(output_file)
        print(f"Saved filtered data to: {output_file}")
        
    def calculate_statistics(self):
        """Calculate statistics with appropriate methods for proteomics"""
        print("\n3. Calculating statistics...")
        
        # Prepare for pairwise comparisons
        self.comparisons = []
        
        # For now, do pairwise comparisons between all groups
        from itertools import combinations
        for group1, group2 in combinations(self.groups, 2):
            print(f"\nComparing {group1} vs {group2}...")
            
            # Get samples for each group
            samples1 = [col for col, g in self.sample_to_group.items() if g == group1]
            samples2 = [col for col, g in self.sample_to_group.items() if g == group2]
            
            # Calculate statistics for each protein
            results = []
            
            for protein_id in self.proteins.index:
                # Get values for each group
                values1 = self.proteins.loc[protein_id, samples1].dropna()
                values2 = self.proteins.loc[protein_id, samples2].dropna()
                
                # Skip if not enough values
                if len(values1) < 2 or len(values2) < 2:
                    continue
                
                # Calculate fold change (on linear scale)
                mean1 = values1.mean()
                mean2 = values2.mean()
                
                # Avoid division by zero or very small numbers
                if mean1 < 1:
                    mean1 = 1
                if mean2 < 1:
                    mean2 = 1
                
                log2fc = np.log2(mean1 / mean2)
                
                # Statistical test - use Welch's t-test (more appropriate for proteomics)
                statistic, pvalue = stats.ttest_ind(values1, values2, equal_var=False)
                
                # Additional statistics
                results.append({
                    'Protein': protein_id,
                    'Gene': self.proteins.loc[protein_id, 'Genes'] if 'Genes' in self.proteins.columns else '',
                    'Protein.Names': self.proteins.loc[protein_id, 'Protein.Names'] if 'Protein.Names' in self.proteins.columns else '',
                    f'Mean_{group1}': mean1,
                    f'Mean_{group2}': mean2,
                    'Log2FC': log2fc,
                    'P.Value': pvalue,
                    f'N_{group1}': len(values1),
                    f'N_{group2}': len(values2)
                })
            
            # Convert to dataframe
            comparison_df = pd.DataFrame(results)
            
            # Multiple testing correction
            if len(comparison_df) > 0:
                _, pvals_corrected, _, _ = multipletests(
                    comparison_df['P.Value'], 
                    alpha=0.05, 
                    method='fdr_bh'
                )
                comparison_df['Q.Value'] = pvals_corrected
                comparison_df['-Log10.P'] = -np.log10(comparison_df['P.Value'])
                
                # Sort by p-value
                comparison_df = comparison_df.sort_values('P.Value')
                
                # Save results
                output_file = self.output_dir / f"comparison_{group1}_vs_{group2}.csv"
                comparison_df.to_csv(output_file, index=False)
                print(f"  Found {(comparison_df['Q.Value'] < 0.05).sum()} significant proteins (Q < 0.05)")
                
                self.comparisons.append({
                    'name': f"{group1}_vs_{group2}",
                    'data': comparison_df,
                    'group1': group1,
                    'group2': group2
                })
                
    def create_visualizations(self):
        """Create various visualizations"""
        print("\n4. Creating visualizations...")
        
        # Set style - use a compatible style
        try:
            plt.style.use('seaborn-darkgrid')
        except:
            plt.style.use('ggplot')  # Fallback style
        sns.set_palette("husl")
        
        # 1. Sample correlation heatmap
        self._plot_correlation_heatmap()
        
        # 2. PCA plot
        if SKLEARN_AVAILABLE:
            self._plot_pca()
        
        # 3. Volcano plots for each comparison
        for comp in self.comparisons:
            self._plot_volcano(comp)
        
        # 4. Heatmap of significant proteins
        self._plot_significant_heatmap()
        
        # 5. Sample quality plots
        self._plot_sample_quality()
        
    def _plot_correlation_heatmap(self):
        """Plot sample correlation heatmap"""
        print("  Creating correlation heatmap...")
        
        # Calculate correlation matrix
        corr_matrix = self.proteins[self.sample_cols].corr()
        
        # Create figure
        plt.figure(figsize=(12, 10))
        
        # Create mask for upper triangle
        mask = np.triu(np.ones_like(corr_matrix), k=1)
        
        # Create heatmap
        sns.heatmap(corr_matrix, 
                   mask=mask,
                   cmap='RdBu_r',
                   center=1,
                   vmin=0.8,
                   vmax=1,
                   square=True,
                   linewidths=0.5,
                   cbar_kws={"shrink": 0.8})
        
        plt.title('Sample Correlation Heatmap', fontsize=16, pad=20)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'correlation_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def _plot_pca(self):
        """Plot PCA"""
        print("  Creating PCA plot...")
        
        # Prepare data
        data_for_pca = self.proteins[self.sample_cols].T.fillna(0)
        
        # Standardize
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data_for_pca)
        
        # PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(data_scaled)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot points colored by group
        for group in self.groups:
            group_samples = [i for i, col in enumerate(self.sample_cols) 
                           if self.sample_to_group.get(col) == group]
            ax.scatter(pca_result[group_samples, 0], 
                      pca_result[group_samples, 1],
                      label=group, s=100, alpha=0.7)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax.set_title('PCA of Samples', fontsize=16)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'pca_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def _plot_volcano(self, comparison):
        """Plot volcano plot for a comparison"""
        print(f"  Creating volcano plot for {comparison['name']}...")
        
        data = comparison['data']
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Define significance thresholds
        fc_threshold = 1  # log2 fold change
        pval_threshold = 0.05  # q-value
        
        # Create color categories
        data['Color'] = 'Not Significant'
        data.loc[(data['Q.Value'] < pval_threshold) & (data['Log2FC'] > fc_threshold), 'Color'] = 'Up'
        data.loc[(data['Q.Value'] < pval_threshold) & (data['Log2FC'] < -fc_threshold), 'Color'] = 'Down'
        
        # Plot
        colors = {'Not Significant': 'gray', 'Up': 'red', 'Down': 'blue'}
        
        for category, color in colors.items():
            mask = data['Color'] == category
            ax.scatter(data.loc[mask, 'Log2FC'], 
                      data.loc[mask, '-Log10.P'],
                      c=color, 
                      alpha=0.6 if category == 'Not Significant' else 0.8,
                      s=20 if category == 'Not Significant' else 40,
                      label=f"{category} ({mask.sum()})")
        
        # Add threshold lines
        ax.axhline(-np.log10(pval_threshold), color='black', linestyle='--', alpha=0.5)
        ax.axvline(fc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(-fc_threshold, color='black', linestyle='--', alpha=0.5)
        
        # Labels
        ax.set_xlabel(f'Log2 Fold Change ({comparison["group1"]}/{comparison["group2"]})')
        ax.set_ylabel('-Log10 Q-Value')
        ax.set_title(f'Volcano Plot: {comparison["name"]}', fontsize=16)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f'volcano_{comparison["name"]}.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save interactive version if plotly available
        if PLOTLY_AVAILABLE:
            self._plot_volcano_interactive(comparison)
            
    def _plot_volcano_interactive(self, comparison):
        """Create interactive volcano plot using plotly"""
        data = comparison['data'].copy()
        
        # Add hover text
        data['hover_text'] = data.apply(
            lambda row: f"Gene: {row['Gene']}<br>" +
                       f"Protein: {row['Protein']}<br>" +
                       f"Log2FC: {row['Log2FC']:.2f}<br>" +
                       f"Q-Value: {row['Q.Value']:.2e}",
            axis=1
        )
        
        # Create figure
        fig = px.scatter(data, 
                        x='Log2FC', 
                        y='-Log10.P',
                        color='Color',
                        hover_name='Gene',
                        hover_data={'hover_text': True, 'Log2FC': False, '-Log10.P': False, 'Color': False},
                        color_discrete_map={'Not Significant': 'gray', 'Up': 'red', 'Down': 'blue'})
        
        # Update layout
        fig.update_layout(
            title=f'Interactive Volcano Plot: {comparison["name"]}',
            xaxis_title=f'Log2 Fold Change ({comparison["group1"]}/{comparison["group2"]})',
            yaxis_title='-Log10 Q-Value',
            hovermode='closest'
        )
        
        # Add threshold lines
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", opacity=0.5)
        fig.add_vline(x=1, line_dash="dash", line_color="black", opacity=0.5)
        fig.add_vline(x=-1, line_dash="dash", line_color="black", opacity=0.5)
        
        # Save
        output_file = self.output_dir / f'volcano_{comparison["name"]}_interactive.html'
        fig.write_html(str(output_file))
        
    def _plot_significant_heatmap(self):
        """Plot heatmap of significant proteins across all comparisons"""
        print("  Creating significant proteins heatmap...")
        
        # Collect all significant proteins
        all_significant = set()
        for comp in self.comparisons:
            significant = comp['data'][comp['data']['Q.Value'] < 0.05]['Protein'].tolist()
            all_significant.update(significant)
        
        if len(all_significant) == 0:
            print("    No significant proteins found for heatmap")
            return
            
        # Subset data to significant proteins
        sig_proteins = self.proteins.loc[list(all_significant), self.sample_cols]
        
        # Z-score normalize
        sig_proteins_z = sig_proteins.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        
        # Create figure
        plt.figure(figsize=(12, max(8, len(sig_proteins) * 0.2)))
        
        # Create clustered heatmap
        sns.clustermap(sig_proteins_z.fillna(0),
                      cmap='RdBu_r',
                      center=0,
                      col_cluster=True,
                      row_cluster=True,
                      cbar_kws={'label': 'Z-score'},
                      figsize=(12, max(8, len(sig_proteins) * 0.2)))
        
        plt.suptitle(f'Heatmap of {len(all_significant)} Significant Proteins', y=0.995)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'significant_proteins_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def _plot_sample_quality(self):
        """Plot sample quality metrics"""
        print("  Creating sample quality plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Number of proteins per sample
        protein_counts = (~self.proteins[self.sample_cols].isna()).sum()
        ax = axes[0, 0]
        protein_counts.plot(kind='bar', ax=ax)
        ax.set_title('Number of Proteins per Sample')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Number of Proteins')
        ax.tick_params(axis='x', rotation=45)
        
        # 2. Intensity distribution per sample
        ax = axes[0, 1]
        self.proteins[self.sample_cols].boxplot(ax=ax, rot=45)
        ax.set_title('Intensity Distribution per Sample')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Log2 Intensity')
        
        # 3. Missing values heatmap
        ax = axes[1, 0]
        missing_data = self.proteins[self.sample_cols].isna().astype(int)
        sns.heatmap(missing_data.T, cmap='RdYlBu', cbar_kws={'label': 'Missing'}, ax=ax)
        ax.set_title('Missing Values Pattern')
        ax.set_xlabel('Proteins')
        ax.set_ylabel('Samples')
        
        # 4. CV distribution by group
        ax = axes[1, 1]
        cv_by_group = {}
        for group in self.groups:
            group_samples = [col for col, g in self.sample_to_group.items() if g == group]
            if len(group_samples) > 1:
                cv_values = self.proteins[group_samples].std(axis=1) / self.proteins[group_samples].mean(axis=1) * 100
                cv_by_group[group] = cv_values.dropna()
        
        pd.DataFrame(cv_by_group).boxplot(ax=ax)
        ax.set_title('CV Distribution by Group')
        ax.set_xlabel('Group')
        ax.set_ylabel('CV (%)')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'sample_quality_plots.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    def prepare_go_analysis(self):
        """Prepare data for GO analysis"""
        print("\n5. Preparing for GO analysis...")
        
        # For each comparison, create gene lists
        for comp in self.comparisons:
            data = comp['data']
            
            # Significant up-regulated genes
            up_genes = data[(data['Q.Value'] < 0.05) & (data['Log2FC'] > 1)]['Gene'].dropna()
            up_genes = [g for g in up_genes if not g.startswith('##')]  # Remove contaminants
            
            # Significant down-regulated genes  
            down_genes = data[(data['Q.Value'] < 0.05) & (data['Log2FC'] < -1)]['Gene'].dropna()
            down_genes = [g for g in down_genes if not g.startswith('##')]
            
            # Save gene lists
            if len(up_genes) > 0:
                output_file = self.output_dir / f"genes_up_{comp['name']}.txt"
                with open(output_file, 'w') as f:
                    f.write('\n'.join(up_genes))
                print(f"  Saved {len(up_genes)} up-regulated genes for {comp['name']}")
                
            if len(down_genes) > 0:
                output_file = self.output_dir / f"genes_down_{comp['name']}.txt"
                with open(output_file, 'w') as f:
                    f.write('\n'.join(down_genes))
                print(f"  Saved {len(down_genes)} down-regulated genes for {comp['name']}")
                
        print("\nGene lists saved. For C. elegans GO analysis, you can use:")
        print("  - WormBase: https://wormbase.org/tools/enrichment")
        print("  - g:Profiler: https://biit.cs.ut.ee/gprofiler/")
        print("  - DAVID: https://david.ncifcrf.gov/")
        
    def create_report(self):
        """Create summary report"""
        print("\n6. Creating summary report...")
        
        report = []
        report.append("DIA-NN Proteomics Analysis Report")
        report.append("=" * 50)
        report.append(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
        report.append(f"\nInput Files:")
        report.append(f"  - Quantification: {self.quant_file}")
        report.append(f"  - Metadata: {self.metadata_file}")
        report.append(f"\nSamples: {len(self.sample_cols)}")
        report.append(f"Groups: {', '.join(self.groups)}")
        report.append(f"Proteins after filtering: {len(self.proteins)}")
        
        report.append(f"\nComparison Results:")
        for comp in self.comparisons:
            data = comp['data']
            report.append(f"\n{comp['name']}:")
            report.append(f"  - Total proteins tested: {len(data)}")
            report.append(f"  - Significant (Q < 0.05): {(data['Q.Value'] < 0.05).sum()}")
            report.append(f"  - Up-regulated (Q < 0.05, FC > 2): {((data['Q.Value'] < 0.05) & (data['Log2FC'] > 1)).sum()}")
            report.append(f"  - Down-regulated (Q < 0.05, FC < 0.5): {((data['Q.Value'] < 0.05) & (data['Log2FC'] < -1)).sum()}")
            
        # Save report
        output_file = self.output_dir / "analysis_report.txt"
        with open(output_file, 'w') as f:
            f.write('\n'.join(report))
            
        print(f"\nAnalysis complete! Results saved to: {self.output_dir}")
        
    def run_analysis(self):
        """Run complete analysis pipeline"""
        self.load_data()
        self.filter_proteins()
        self.calculate_statistics()
        self.create_visualizations()
        self.prepare_go_analysis()
        self.create_report()

def main():
    if len(sys.argv) < 3:
        print("Usage: python diann_analysis.py <maxlfq_file> <metadata_file>")
        print("\nExample:")
        print("  python diann_analysis.py DR12_MaxLFQ.csv DR12_metadata.csv")
        sys.exit(1)
        
    quant_file = sys.argv[1]
    metadata_file = sys.argv[2]
    
    # Run analysis
    analyzer = DIANNAnalyzer(quant_file, metadata_file)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()