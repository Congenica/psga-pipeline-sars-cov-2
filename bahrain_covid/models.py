# coding: utf-8
from sqlalchemy import Boolean, Column, DateTime, Enum, Float, ForeignKey, Integer, LargeBinary, String, Table, Text, text
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class Area(Base):
    __tablename__ = 'area'
    __table_args__ = {'schema': 'sars_cov_2', 'comment': 'Geographic "Area" location'}

    name = Column(String, primary_key=True, comment='Primary key, the area name, capitol level geolocation from where the sample was collected')
    latitude = Column(Float(53), comment='The angular distance of the location north or south of the equator')
    longitude = Column(Float(53), comment='the angular distance of the location east or west of the Greenwich meridian')
    colour = Column(LargeBinary, comment='The colour of this block in reports')


class Block(Base):
    __tablename__ = 'block'
    __table_args__ = {'schema': 'sars_cov_2', 'comment': 'Geographic "Block" location'}

    number = Column(String, primary_key=True, comment='Primary key, the block number, low level geolocation from where the sample was collected')
    latitude = Column(Float(53), comment='The angular distance of the location north or south of the equator')
    longitude = Column(Float(53), comment='the angular distance of the location east or west of the Greenwich meridian')
    colour = Column(LargeBinary, comment='The colour of this block in reports')


class Comorbidity(Base):
    __tablename__ = 'comorbidity'
    __table_args__ = {'schema': 'sars_cov_2', 'comment': 'Coded diseases presented in patient (eg Diabetes)'}

    comorbidity_id = Column(Integer, primary_key=True, comment='Primary key')
    description = Column(Text, comment='Description of the comorbidity')

    samples = relationship('Sample', secondary='sars_cov_2.sample_comorbidity')


class Governorate(Base):
    __tablename__ = 'governorate'
    __table_args__ = {'schema': 'sars_cov_2', 'comment': 'Geographic "Governorate" location'}

    name = Column(Enum('Capital', 'Northern', 'Southern', 'Muharraq', name='governorate_name'), primary_key=True)
    latitude = Column(Float(53))
    longitude = Column(Float(53))
    colour = Column(LargeBinary)


class Sample(Base):
    __tablename__ = 'sample'
    __table_args__ = {'schema': 'sars_cov_2', 'comment': 'Sample data'}

    sample_id = Column(Integer, primary_key=True, server_default=text("nextval('\"sars_cov_2\".sample_sample_id_seq'::regclass)"), comment='Primary key, serial sequence generated in the database')
    lab_id = Column(String, comment='Lab sample identifier')
    date_collected = Column(DateTime, comment='timestamp of sample collection')
    data_sequenced = Column(DateTime, comment='timestamp of sample sequencing')
    sequencing_run = Column(String, comment='identifier of sequencing run')
    gender = Column(Enum('M', 'F', 'U', name='gender'), nullable=False, server_default=text("'U'::gender"), comment='Gender of host person')
    age = Column(Integer, comment='Age of host person in years')
    nationality = Column(String, comment='Nationality of host person')
    governorate_name = Column(ForeignKey('sars_cov_2.governorate.name', deferrable=True, initially='DEFERRED'), comment='high level geolocation from where sample was collected')
    area_name = Column(ForeignKey('sars_cov_2.area.name', deferrable=True, initially='DEFERRED'), comment='foreign key to area table, mid level geolocation from where sample was collected')
    block_number = Column(ForeignKey('sars_cov_2.block.number', deferrable=True, initially='DEFERRED'), comment='foreign key to block table, low level geolocation from where sample was collected')
    sample_number = Column(Integer, comment='Lab sample identifier (may be redundant with lab_id)')
    ct_value = Column(Float(53), comment='RNA quality measure')
    symptoms = Column(Text, comment='Symptoms observed')
    travel_exposure = Column(String, comment='Country host visited recently')
    hospital_admittance = Column(Enum('No', 'Standard', 'ICU', name='hospital_admittance'), comment='Hospital type host admited to')
    pangolin_lineage = Column(String, comment='Viral phylogenetic lineage')
    genbank_id = Column(Integer, comment='GeneBank accession of virus genome')
    gisaid_id = Column(String, comment='GISAID identifier (GISAID is a global science initiative and primary source that provides open-access to genomic data of the novel coronavirus responsible for COVID-19)')
    genome_length = Column(Integer, comment='Number of base pair in the virus genome')
    mrn = Column(Integer, comment='Individual National Identity number')

    area = relationship('Area')
    block = relationship('Block')
    governorate = relationship('Governorate')


t_sample_comorbidity = Table(
    'sample_comorbidity', metadata,
    Column('sample_id', ForeignKey('sars_cov_2.sample.sample_id', deferrable=True, initially='DEFERRED'), primary_key=True, nullable=False, comment='Primary key with comorbidity_id, and foreign key to sample table'),
    Column('comorbidity_id', ForeignKey('sars_cov_2.comorbidity.comorbidity_id', deferrable=True, initially='DEFERRED'), primary_key=True, nullable=False, comment='Primary key with sample_id, and foreign key to comorbidity table'),
    schema='sars_cov_2',
    comment='Sample comorbitiy, a linking table between the sample an comorbidity record'
)


t_sample_qc = Table(
    'sample_qc', metadata,
    Column('sample_id', ForeignKey('sars_cov_2.sample.sample_id', deferrable=True, initially='DEFERRED'), comment='Primary key, and foreign key to sample table'),
    Column('pct_N_bases', Float(53), comment='percentage of N bases in sequenced genome'),
    Column('pct_covered_bases', Float(53), comment='percentage of covered bases in sequenced genome'),
    Column('longest_no_N_run', Integer, comment='longest sequence run without N base'),
    Column('num_aligned_reads', Float(53), comment='Number of reads succesfully aligned and mapped'),
    Column('qc_pass', Boolean, comment='Has sample passed QC'),
    Column('qc_plot', LargeBinary, comment='QC plot image'),
    Column('pipeline_version', String, comment='mapping pipeline version'),
    schema='sars_cov_2',
    comment='Sample Quality Control, QC data for a sample'
)
