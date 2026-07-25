import csv, json
from pathlib import Path
base = Path(r'E:\gradthesis\curvecurator\8171\GSEA_from_raw_name_pEC50')
rank_paths = {
    'best_pEC50': base/'rank_best_pEC50.csv',
    'weighted_score': base/'rank_weighted_score.csv',
    'mean_pEC50': base/'rank_mean_pEC50.csv',
}
map_path = base/'symbol_to_entrez_mapping.csv'
driver_paths = [
    base/'weighted_score'/'selected_pathway_driver_genes.csv',
    base/'weighted_score'/'top5_pathway_driver_genes.csv',
    base/'weighted_score'/'top5_pathway_driver_genes_with_scores.csv',
]
# mapping
entrez_to_symbol = {}
symbol_to_entrez = {}
with open(map_path, newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        e = row.get('ENTREZID') or row.get('entrez') or row.get('ENTREZ')
        s = row.get('SYMBOL') or row.get('symbol')
        if e and s:
            entrez_to_symbol[str(e)] = s
            symbol_to_entrez[s] = str(e)
# ranks
ranks = {}
for name, path in rank_paths.items():
    rows = []
    with open(path, newline='', encoding='utf-8-sig') as f:
        for i, row in enumerate(csv.DictReader(f), start=1):
            e = str(row.get('ENTREZID') or row.get('ENTREZ') or row.get('gene') or '')
            score = row.get('score') or row.get('stat') or row.get(name) or ''
            rows.append((e, i, score, entrez_to_symbol.get(e, e)))
    ranks[name] = {e: {'rank': i, 'score': score, 'symbol': sym} for e, i, score, sym in rows}
    ranks[name]['__n__'] = len(rows)
# driver symbols/entrez
interesting = ['NUP205','TPR','DHX9','RANBP2','XPO1','IPO5','KPNA2','NUP98','TNPO1','LRPPRC','CKAP5']
# add from driver files
for path in driver_paths:
    if path.exists():
        with open(path, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                s = row.get('SYMBOL') or row.get('symbol') or row.get('gene') or row.get('Gene')
                e = row.get('ENTREZID') or row.get('entrez') or row.get('ENTREZ')
                if s and s not in interesting: interesting.append(s)
                if e and e in entrez_to_symbol and entrez_to_symbol[e] not in interesting: interesting.append(entrez_to_symbol[e])
# summarize unique first 20 with ranks in best/weighted
summary = []
seen=set()
for s in interesting:
    if s in seen: continue
    seen.add(s)
    e = symbol_to_entrez.get(s)
    if not e: continue
    item = {'symbol': s, 'entrez': e}
    for name in ['best_pEC50','weighted_score','mean_pEC50']:
        r = ranks[name].get(e)
        if r:
            item[name+'_rank']=r['rank']
            item[name+'_score']=r['score']
            item[name+'_pct']=round(r['rank']/ranks[name]['__n__']*100,1)
    summary.append(item)
# Sort by weighted rank if available else best
summary.sort(key=lambda x: x.get('weighted_score_rank', x.get('best_pEC50_rank', 999999)))
out = Path(r'C:\Users\owner\Documents\Codex\2026-07-07\ni-ha\outputs\real_gsea_driver_ranks.json')
out.write_text(json.dumps({'n_best':ranks['best_pEC50']['__n__'],'n_weighted':ranks['weighted_score']['__n__'],'drivers':summary[:30]}, indent=2), encoding='utf-8')
print(out)
print(json.dumps({'n_best':ranks['best_pEC50']['__n__'],'n_weighted':ranks['weighted_score']['__n__'],'drivers':summary[:15]}, indent=2))
