function M1 = M1_matrix(Db)

    Dbx = Db(1);
    Dby = Db(2);
    Dbz = Db(3);

    M1 = [   Dbx,     0,     Dbz,    - Dby;
               0,   Dbx,     Dby,      Dbz;
             Dbz,   Dby,   - Dbx,        0;
           - Dby,   Dbz,       0,    - Dbx];

end