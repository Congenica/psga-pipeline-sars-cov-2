Submit-block ::= {
  contact {
    contact {
      name name {
        last "LastName",
        first "FirstName",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "MyCompanyName",
        div "MyDivision",
        city "MyCity",
        country "MyCountry",
        street "MyStreet",
        email "my_email@domain.br",
        postal-code "ZIP CODE"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "LastName",
            first "FirstName",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "MyCompanyName",
        div "MyDivision",
        city "MyCity",
        country "MyCountry",
        street "MyStreet",
        postal-code "ZIP CODE"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "LastName",
              first "FirstName",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "Test publication"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:my_email@domain.br"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
