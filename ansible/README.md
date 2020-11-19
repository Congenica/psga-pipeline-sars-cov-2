# Congenica SARS-CoV-2 pipeline install

The anible playbook 'pipeline_install.yml' will install and configure an environment as defined in [SAP-18290](https://jira.congenica.net/browse/SAP-18290).

The following variables need to be set for the playbook to run:

- `_pipeline_user` : the user that will be configured to run the pipeline, i.e 'pipeline'
- `_github_repo` : the name of the GitHub repo to be cloned, i.e 'Congenica/ps-bahrain-covid'
- `_github_checkout_dir` : the directory that the github repo should be cloned into, i.e '~/ps-bahrain-covid' (the default home (~) directory is for the `_pipeline_user` user)
- `_s3_bucket` : the S3 bucket that should be synchronised, i.e 's3://congenica-development-data-share/SAP-18211_Bahrain_COVID'
- `_s3_sync_dir` : the directory used for the S3 sync, i.e '~/SAP-18211_Bahrain_COVID' (the default home (~) directory is for the `_pipeline_user` user)

These are configured inside the 'pipeline_install.yml' file istelf.

When running the playbook, the user will be prompted to set some additional variables:

- `_github_username` : username of the account that will be used to clone the specified repo
- `_github_password_or_token` : password/token of the account that will be used to clone the specified repo
- `_aws_key_id` : This keypair will need permission to sync the specified S3 bucket
- `_aws_private_key` : This keypair will need permission to sync the specified S3 bucket

The playbook can be run against the local machine with the following command:

```shell
ansible-playbook --connection=local -i 127.0.0.1, ansible/pipeline_install.yml
```

variables can also be passed through the command. These will override anything set in the playbook...

```shell
ansible-playbook --connection=local -i 127.0.0.1, ansible/pipeline_install.yml --extra-vars "_pipeline_user=foo"
```

Environment variables can be set for the `_pipeline_user` user by adding them to `vars/env_vars.yml`...

```shell
environment_vars:
  - key: GENOME_FASTA_PATH
    value : '~/foo'
  - key: DIFFERENT_VAR
    value : '~/bar'
```
