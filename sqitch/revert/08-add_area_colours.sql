-- Revert ps-bahrain-covid:08-add_area_colours from pg

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  UPDATE area SET colour = NULL WHERE name IN (
    'ZINJ',
    'ZAYED TOWN',
    'ZALLAQ',
    'WADIYAN',
    'WADI ALI',
    'UMMHAZWARAH',
    'UMM JIDR AL SUMMAN',
    'UMM JIDR',
    e'UMM AL NA\'SAN',
    'UM AL BAIDH',
    'TUBLI',
    'TASHAN',
    'TARRAFI',
    'SUWAD AL SHAMALIYAH',
    'SUWAD AL JANUBIYAH',
    'SUFALA',
    'SOUTH SEHLA',
    'SHARQ AL HIDD UM AL SHAJAR',
    'SHARQ AL HIDD UM AL SAWALI',
    'SHARQ AL HIDD AL QULYAH',
    'SHARQ AL HIDD AL EZEL',
    'SHAKHURAH',
    'SHAHRAKKAN',
    'SAR',
    'SANAD',
    'SANABIS',
    'SAMAHEEJ',
    'SALMAN TOWN',
    'SALIMABAD',
    'SALIHIYA',
    'SAFREH',
    'SADAD',
    'RUBUD AL SHARQIYAH',
    'RUBUD AL GHARBIYAH',
    'RIFFA MOASKAR',
    'RIFFA / WADI ALSAIL',
    'RIFFA / SWAYFRA',
    e'RIFFA / MO\'ASKAR',
    'RIFFA / JARYALSHAIKH',
    'RIFFA / BUKOWARAH',
    'RIFFA / ALSHARGI',
    'RIFFA / ALROWDHA',
    'RIFFA / ALHUNAYNIYAH',
    'RIFFA / ALHAJIYAT',
    'RIFFA / ALGHARBI',
    'RIFFA / AL SHAMALI',
    'RIFFA / AL BUHAIR',
    'RAS ZUWAYED',
    'RAS HAYYAN',
    'RAS ABU JARJUR',
    'QALALI',
    'NURANA ISLAND',
    'NORTH SEHLA',
    'NABIH SALEH',
    'MUZARRA',
    'MURQOBAN',
    'MUHARRAQ',
    'MOHAMMADIYAH',
    'MAQABAH',
    'MANAMA/MINA SALMAN INDUSTRIAL AREA',
    'MANAMA CENTER',
    'MANAMA / UMM ALHASSAM',
    'MANAMA / SEA FRONT',
    'MANAMA / BUGHAZAL',
    'MANAMA / BUASHIRAH',
    'MANAMA / ALSUWAYFIYAH',
    'MANAMA / ALSUQAYYAH',
    'MANAMA / ALSALMANIYA',
    'MANAMA / ALQUDAYBIYAH',
    e'MANAMA / ALNAI\'M',
    'MANAMA / ALMAHUZ',
    'MANAMA / ALJUFFAIR',
    'MANAMA / ALHOORA',
    'MANAMA / ALGUFUL',
    'MANAMA / ALGHURAYFAH',
    'MANAMA / ALFATEH',
    'MANAMA / ALCORNICHE',
    'MANAMA / AL ADLIYAH',
    'MAMLAHAT AL MAMTALAH',
    'MALKIYA',
    'MAHAZZAH',
    'MADINAT KHALIFA ALULAIM',
    'MADINAT KHALIFA ALMAHADIR',
    'MADINAT HAMAD',
    'LHASSAY',
    'KING FAHAD CAUSWAY',
    'KHAMIS',
    'KARZAKKAN',
    'KARRANAH',
    'KARBABAD',
    'JUZUR AL DAR',
    'JIDHAFS',
    'JIDDAH',
    'JIDD ALI',
    'JID AL HAJ',
    'JERDAB',
    'JEBLAT HIBSHI',
    'JAZIRAT HAWAR',
    'JAZAER BEACH',
    'JAU',
    'JANNUSAN',
    'ISA TOWN',
    'INDUSTRIAL AREA',
    'HILLAT ABDULSALEH',
    'HIDD AL JAMAL',
    'HIDD',
    'HAWRAT SANAD',
    'HAWRAT INGAH',
    'HAWRAT A ALI',
    'HALAT AL SLETAH',
    'HALAT AL NEAIM',
    'HAFIRAH',
    'DURRAT AL BAHRAIN',
    'DIYAR AL MUHARRAQ',
    'DIPLOMATIC AREA',
    'DILMUNIA ISLAND',
    'DAR KULAIB',
    'DAMISTAN',
    'BUSAITEEN',
    'BURI',
    'BUDAIYA',
    'BU QUWAH',
    'BILAD AL QADEEM',
    'BAR BAR',
    'BANI JAMRAH',
    'BANDER AL SEEF',
    'Area',
    'AWALI',
    'ASKAR',
    'ARAD',
    'AMWAJ',
    'AL SHABAK',
    'AL SEEF',
    'AL SAYH',
    'AL SAKHIR',
    'AL SAFRIYAH',
    'AL RUMAYTHAH',
    'AL RUMAMIN',
    'AL RIFFAH',
    'AL RAMLI',
    'AL QURAYYAH',
    'AL QURAYN',
    'AL QARYAH',
    'AL QARAH',
    e'AL QAL\'AH',
    'AL QADAM',
    'AL NUWAIDRAT',
    'AL NASFA',
    'AL MUSALLA',
    'AL MAZROWIAH',
    'AL MARKH',
    'AL MAQSHA',
    'AL MAMTALAH',
    e'AL MA\'AMEER',
    'AL LAWZI',
    'AL KHARIJIYAH',
    'AL JASRAH',
    'AL JANABIYAH',
    e'AL JA\'SIRAH',
    'AL HAMRIYA',
    'AL HAMALAH',
    'AL HAJAR',
    'AL GHAYNAH',
    'AL DUR',
    'AL DIRAZ',
    'AL DAIR',
    'AL DAIH',
    'AL BURHAMA',
    'AL AMUR',
    'AL AKR AL SHARQI',
    'AL AKR AL GHARBI',
    'ADHARI',
    'ABU SAYBA',
    'ABU BHAM',
    e'ABU AL \'AYASH',
    e'A\'ALI'
  );

COMMIT;
