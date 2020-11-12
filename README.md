# ps-bahrain-covid

## Sqitch

`Sqitch` manages all schema migrations. Prerequisites for running the `sqitch`:

* `sqitch` is installed in the machine. See [sqitch downloads page](https://sqitch.org/download/)
* `Postgres` with `psql` installed in the environment. Has postgres user set up
* `Postgres` has password exported to env variable 
```commandline
export PGPASSWORD=some_secret_password
```
* A dedicated database is created for the project. It may be created using the following cmd:
```commandline
createdb -h localhost -U postgres bahrain_sars_cov_2
```

### Working with sqitch

The work is done in `sqitch/` directory

To add the new migration, use the following command:
```commandline
sqitch add 0x-my_migration_file_name -n 'My changes described here'
```
An `.sql` file will be created in `deploy`, `revert` and `verify`. Populate the files.
Important - `verify` scripts only fail, if `.sql` query raises an exception. Create verify scripts
accordingly to throw exceptions in case of failed verification

To check the migration status (are we missing any migrations?):
```commandline
sqitch status db:pg:bahrain_sars_cov_2
```

Migrations can be made using the following:
```commandline
sqitch deploy db:pg:bahrain_sars_cov_2
```

To verify migrations, which were made:
```commandline
sqitch verify db:pg:bahrain_sars_cov_2
```

To revert the changes:
```commandline
sqitch revert db:pg:bahrain_sars_cov_2
```

### Sqitch troubleshoot

If authentication fails, try adding connection info to the command-line. For example:
```commandline
sqitch --db-user postgres --db-host localhost --db-port 5432 deploy db:pg:bahrain_sars_cov_2
```