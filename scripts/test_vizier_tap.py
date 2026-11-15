import pyvo

try:
    tap = pyvo.dal.TAPService('https://tapvizier.cds.unistra.fr/TAPVizieR/tap')
    q = 'SELECT TOP 3 * FROM "II/328/allwise"'
    r = tap.search(q)
    print('rows', len(r))
    print(r.colnames[:10])
except Exception as e:
    print('ERR', str(e)[:300])
