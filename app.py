from flask import Flask, request, url_for, redirect, render_template, send_from_directory
from ml_models import get_prediction
from seq_sim import *
from struct_sim import *
from flask_wtf import FlaskForm
from wtforms import FileField, TextAreaField, StringField, IntegerField, validators, SubmitField
from Bio import SeqIO
import subprocess

from flask_bootstrap import Bootstrap5


app = Flask(__name__)
app.config['SECRET_KEY'] = 'you-shall-not-pass'
bootstrap = Bootstrap5(app)

foldseek_initcode = ['tar xvzf assets/foldseek-linux-avx2.tar.gz'.split(), 'export PATH=$(pwd)/foldseek/bin/:$PATH'.split()]
for initcode in foldseek_initcode:
    print(initcode)
    subprocess.run(initcode, shell=True)

# Custom validator for fasta_string and fasta_file

class UploadForm(FlaskForm):
    fasta_string = TextAreaField("Input sequence", render_kw = {"placeholder": "Paste your FASTA sequence here"}, validators=[validators.optional()])
    fasta_file = FileField("Or upload FASTA file", validators=[validators.optional()])
    email_return = StringField("(Optional) E-mail to return results", validators=[validators.email(), validators.optional()])
    avp_threshold = IntegerField("AVP prediction threshold", default = 50, validators=[validators.input_required(), validators.NumberRange(min=0, max=100, message='Validators must be within 0-100')])
    anticov_threshold = IntegerField("ACoVP prediction threshold", default = 50, validators=[validators.input_required(), validators.NumberRange(min=0, max=100, message='Validators must be within 0-100')])
    submit_button = SubmitField()

    def validate(self):
        if not super(UploadForm, self).validate():
            return False
        if not(bool(self.fasta_string.data) ^ bool(self.fasta_file.data)):
            self.fasta_string.errors += (validators.ValidationError('Please either enter sequence or upload file'),)
            return False
        return True

def process(option_dict: dict = None) -> None :
    all_fasta = {'id':[], 'seq':[]}
    with open('input/uploaded.fasta') as fhandle:
        for record in SeqIO.parse(fhandle, 'fasta'):
            all_fasta['id'].append(str(record.id))
            all_fasta['seq'].append(str(record.seq))
    # Run ML Model
    get_prediction(all_fasta, avp_threshold = option_dict['avp_threshold'], anticov_threshold = option_dict['anticov_threshold'], rename_dict = {'SeqID': 'Sequence ID', 'prediction_label': 'Predicted Class', 'prediction_score': 'Probability'})
    
    # Run sequence sim model #TODO: Move into seq_sim.py
    get_blast_result("uploaded.fasta")

    # Run structure sim model #TODO: Move into struct_sim.py
    structure_process(all_fasta)

def generate_html_features(result_path: str, col_count: int) -> None:
    data = pd.read_csv(result_path).loc[:, :col_count]
    data_columns = data.columns    
    data_values = data.values.tolist()
    return data_columns, data_values

@app.route('/', methods=['GET', 'POST'])
def upload_func():
    form = UploadForm()
    if request.method == 'POST' and form.validate():
        if bool(form.fasta_file.data):
            form.fasta_file.data.save('input/uploaded.fasta')
        if bool(form.fasta_string.data):
            with open('input/uploaded.fasta', 'w') as inputhandle:
                inputhandle.write(form.fasta_string.data)
        option_dict = {'email_return': form.email_return.data, #TODO: Fix if user does not specify email
                       'avp_threshold': form.avp_threshold.data/100,
                       'anticov_threshold': form.anticov_threshold.data/100}
        process(option_dict)
        return redirect(url_for('avp_result_func'))
    return render_template('upload.html', form = form)

@app.route('/about')
def about_func():
    return render_template('about.html')

@app.route('/avp_result')
def avp_result_func():
    data_columns, data_values = generate_html_features('output/avp_result.csv', col_count = 3)
    # A dictionary of feature explanation {feature_name: feature_explanation}
    feature_explanation = {'Probability': 'Probability of being an AVP'}
    return render_template('result.html', 
                           page_header = "AVP: Machine learning prediction result",
                           header_list = data_columns, 
                           data = data_values, 
                           explanation = feature_explanation)

@app.route('/anticov_result')
def anticov_result_func():
    data_columns, data_values = generate_html_features('output/anticov_result.csv', col_count = 3)
    return render_template('result.html', 
                           page_header = "ACoVP: Machine learning prediction result",
                           header_list = data_columns, 
                           data = data_values, 
                           explanation = {})

@app.route('/avp_blast_result')
def avp_blast_result_func():
    data_columns, data_values = generate_html_features('output/avp_blast_output.csv', col_count = 5)
    feature_explanation = {'Class': 'Peptide class according to the best matched sequence from our database (AVP/Non-AVP)'}
    return render_template('result.html', 
                           page_header = "AVP: Sequence similarity result",
                           header_list = data_columns, 
                           data = data_values, 
                           explanation = feature_explanation)

@app.route('/anticov_blast_result')
def anticov_blast_result_func():
    data_columns, data_values = generate_html_features('output/anticov_blast_output.csv', col_count = 5)
    feature_explanation = {'Class': 'Peptide class according to the best matched sequence from our database (Anti-CoV/AVP)'}
    return render_template('result.html', 
                           page_header = "ACoVP: Sequence similarity result",
                           header_list = data_columns, 
                           data = data_values, 
                           explanation = feature_explanation)

@app.route('/avp_structure_result')
def avp_structure_result_func():
    data_columns, data_values = generate_html_features('output/avp_structure_result.csv', col_count = 5)
    feature_explanation = {'Class': 'Peptide class according to the best matched sequence from our database (Anti-CoV/AVP)'}
    return render_template('result.html',
                           page_header = "AVP: Structure similarity result",
                           header_list = data_columns, 
                           data = data_values, 
                           explanation = feature_explanation)

@app.route('/anticov_structure_result')
def anticov_structure_result_func():
    data_columns, data_values = generate_html_features('output/anticov_structure_result.csv', col_count = 5)
    feature_explanation = {'Class': 'Peptide class according to the best matched sequence from our database (Anti-CoV/AVP)'}
    return render_template('result.html',
                           page_header = "ACoVP: Structure similarity result",
                           header_list = data_columns, 
                           data = data_values, 
                           explanation = feature_explanation)

@app.route('/download_avp_result')
def download_avp_result_func():
    return send_from_directory('output', 'avp_result.csv', as_attachment=True) 

@app.route('/download_anticov_result')
def download_anticov_result_func():
    return send_from_directory('output', 'anticov_result.csv', as_attachment=True) 


    