function M2 = M2_matrix(Db)

    Dbx = Db(1);
    Dby = Db(2);
    Dbz = Db(3);

    M2 = [   Dby, - Dbz,    0,    Dbx;
           - Dbz, - Dby,  Dbx,      0;
               0,   Dbx,  Dby,    Dbz;
             Dbx,     0,  Dbz,  - Dby];
end