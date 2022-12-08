function M3 = M3_matrix(Db)

    Dbx = Db(1);
    Dby = Db(2);
    Dbz = Db(3);

    M3 = [   Dbz,   Dby, - Dbx,    0;
             Dby, - Dbz,     0,  Dbx;
           - Dbx,     0, - Dbz,  Dby;
               0,   Dbx,   Dby,  Dbz];
end