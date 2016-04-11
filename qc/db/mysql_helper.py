def last_insert_id(mysql_cursor):
    mysql_cursor.execute('SELECT LAST_INSERT_ID() AS id')
    m_row = mysql_cursor.fetchone()
    if m_row:
        return m_row[0]
    raise Exception("No id found from LAST_INSERT_ID()")
