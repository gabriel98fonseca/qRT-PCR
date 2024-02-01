import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import streamlit as st
import io

st.set_page_config(page_title="RT-qPCR Analysis")
st.title('RT-qPCR Analysis')
st.header("Upload your CSV file")
file_uploaded = st.file_uploader("Upload a CSV file", type=["csv", "xlsx"], label_visibility="hidden")

if file_uploaded:
    df = pd.read_csv(file_uploaded, skiprows=6, header=1, usecols=[1, 2, 9], names=['Sample', 'Target', 'Ct'])
    df['Ct'] = df['Ct'].astype(float)

    # Sort the DataFrame by "Target" and "Sample"
    df = df.sort_values(by=['Target', 'Sample']).reset_index(drop=True)

    # Display the DataFrame with checkboxes using st.experimental_data_editor
    st.header("Remove outliers")
    edited_df = st.data_editor(df, num_rows='dynamic', width=750)

    st.header("Check your Standard Deviation")
    merge = edited_df.groupby(['Sample', 'Target'])['Ct'].agg(['mean', 'std']).reset_index()

    # Display checkboxes for each sample
    selected_sample_list = st.sidebar.multiselect("Select Samples:", merge['Sample'].unique())

    if selected_sample_list:
        # Filter the dataframe based on selected samples
        merge = merge[merge['Sample'].isin(selected_sample_list)]

        st.dataframe(merge)
        with st.sidebar:
            Samp = st.selectbox("Select your reference sample:", merge['Sample'].unique())
            reference_gene = st.selectbox('Select your reference gene:', merge['Target'].unique())
            Tar = merge['Target'].unique()


        def calculate_ddCt(df, targets, reference_gene, control_sample):
            final_ddCt_data = pd.DataFrame()

            for target_gene in targets:
                # Skip processing for the reference gene
                if target_gene == reference_gene:
                    continue

                target_data = df[df['Target'] == target_gene]
                reference_data = df[df['Target'] == reference_gene]

                merged_data = pd.merge(target_data, reference_data, on='Sample', suffixes=('_target', '_reference'))

                # Calculate dCt (difference in Ct values)
                merged_data['dCt'] = merged_data['mean_target'] - merged_data['mean_reference']

                # Calculate ddCt (difference in dCt values relative to a control condition)
                control_dCt = merged_data[merged_data['Sample'] == control_sample]['dCt'].values[0]
                merged_data['ddCt'] = merged_data['dCt'] - control_dCt

                # Normalizing data
                merged_data['norm'] = 2 ** -merged_data['ddCt']

                # Calculate combined error

                merged_data['combined error'] = ((merged_data['std_reference']**2)+(merged_data['std_target']**2))**0.5

                # Calculating RQ values

                merged_data['RQ min'] = merged_data['norm'] - 2** -(merged_data['ddCt']+ merged_data['combined error'])

                merged_data['RQ max'] = 2**-(merged_data['ddCt'] - merged_data['combined error']) - merged_data['norm'] 

                # Concatenate results to the final DataFrame
                final_ddCt_data = pd.concat([final_ddCt_data, merged_data[['Sample', 'Target_target', 'norm','RQ min', 'RQ max']]]).reset_index(drop=True)

            return final_ddCt_data

        final_df = calculate_ddCt(merge, Tar, reference_gene, Samp)

        # Plot a bar graph for each target with samples as columns and RQ min/max as errors
            
        # Dropdown to select target gene
        selected_target_gene = st.selectbox("Select Target Gene", Tar)

        # Check if a file has been uploaded
        if selected_target_gene is not None:
            if selected_target_gene == reference_gene:
                st.warning("Please select a non-reference target gene.")
            else:
                # Filter data for the selected target gene
                target_data = final_df[final_df['Target_target'] == selected_target_gene]

                # Plot the selected target gene
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.bar(target_data['Sample'], target_data['norm'], yerr=[target_data['RQ min'], target_data['RQ max']], capsize=5, label=selected_target_gene)
                ax.set_xlabel('Sample')
                ax.set_ylabel('Normalized Expression')
                ax.set_title(f'Expression of {selected_target_gene}')

                # Set the background of the plot to be transparent
                ax.patch.set_alpha(0.0)

                # Remove outer border
                ax.spines['top'].set_color('none')
                ax.spines['right'].set_color('none')

                # Display the plot in the Streamlit app
                st.pyplot(fig)

                # Button to download the plot as a PDF
            download_button = st.button("Download PDF")

            if download_button and plt.fignum_exists(fig.number):
            # Save the plot to a BytesIO buffer
                buffer = io.BytesIO()
                plt.savefig(buffer, format='pdf', bbox_inches='tight', transparent=True)
                buffer.seek(0)

                # Create a download link for the PDF
                st.download_button(label="Download PDF", data=buffer, file_name=f"{selected_target_gene}.pdf", key="download_pdf")
