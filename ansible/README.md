# Congenica SARS-CoV-2 pipeline install

The anible playbook 'pipeline_install.yml' will install and configure an environment as defined in [SAP-18290](https://jira.congenica.net/browse/SAP-18290).

### Configure Ansible
As the user who will perform the install ensure that ansible is installed and functional

The following variables need to be set for the playbook to run. The user will be prompted to enter these when executing the ansible playbook:

 - `db_host` : The host of the database, default: localhost
 - `db_name` : The logical name of the database on the host, default: bahrain_sars_cov_2
 - `db_user` : The user to use to connect to the database, default: postgres
 - `db_port` : The port to use to connect to the database, default: 5432
 - `db_password` : The password to use to connect to the database. No default.

 - `covid_pipeline_rootdir` :Path to the pipeline code (e.g. git checkout). Default: ${HOME}/ps-bahrain-covid
 - `covid_pipeline_fastq_path` : Path to the input FASTQ files and TSV metadata file. Default: ${HOME}/Bahrain_COVID_s3_data_lite/sample_data
 - `covid_pipeline_workdir` :Path to the whole pipeline output. Default: ${HOME}/covid-pipeline
 - `covid_pipeline_reports_path` :   Path to the pipeline reports. Default: ${HOME}/reports


The playbook can be run against the local machine with the following command:


```shell
ansible-playbook --connection=local -i 127.0.0.1, ansible/pipeline_install.yml
```

variables can also be passed through the command. These will override anything set in the playbook...

```shell
ansible-playbook --connection=local -i 127.0.0.1, ansible/pipeline_install.yml --extra-vars "_pipeline_user=foo"
```

Environment variables can be set for the installing user by adding them to `vars/env_vars.yml`...

> Ensure that any additiona vars and values are set here before running the playbook.

```shell
environment_vars:
  - key: GENOME_FASTA_PATH
    value : '~/foo'
  - key: DIFFERENT_VAR
    value : '~/bar'
```

## Troubleshooting
### ansible-galaxy or ansible-playbook are not installed.
Ensure ansible is installed with the commands below, run from the user attempting the install:
```shell
sudo yum install python3 python3-pip -y
pip3 install ansible --user --upgrade
```
### ansible-playbook fails with error for missing resource docker_image
You may receive an error:
```shell
ERROR! couldnt resolve module/action community.general.docker_image. This often indicates a misspelling, missing collection, or incorrect module path.
```
Manually install dependencies:
```shell
ansible-galaxy collection install -r ansible/requirements.yml
```

### Post install failures to communicate with docker
You may receive an error such as:
```shell
Got permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get http://%2Fvar%2Frun%2Fdocker.sock/v1.24/containers/json: dial unix /var/run/docker.sock: connect: permission denied
```

The ansible playbook modifies the running users groups, but this will not take effect until the next new session. Log out of your current session and log back in. Running `groups` should now show a secondary group of `docker`