import sys
import os
import json
import logging
from bioagents import Bioagent
from .causality_agent import CausalityAgent
from indra.sources.trips.processor import TripsProcessor
from kqml import KQMLModule, KQMLPerformative, KQMLList, KQMLString, KQMLToken
from indra.statements import stmts_from_json, Agent
from indra.databases import hgnc_client


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('CausalA')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/resources/'


class CausalityModule(Bioagent):
    name = 'CausalA'
    tasks = ['FIND-CAUSAL-PATH', 'FIND-CAUSALITY-TARGET',
             'FIND-CAUSALITY-SOURCE',
             'DATASET-CORRELATED-ENTITY', 'FIND-COMMON-UPSTREAMS',
             'RESTART-CAUSALITY-INDICES', 'FIND-MUTEX', 'FIND-MUTATION-SIGNIFICANCE',
             'RESET-CAUSALITY-INDICES',  'FIND-CELLULAR-LOCATION-FROM-NAMES',
             'FIND-CELLULAR-LOCATION', 'FIND-GENE-SUMMARY']

    def __init__(self, **kwargs):
        self.CA = CausalityAgent(_resource_dir)
        # Call the constructor of KQMLModule
        super(CausalityModule, self).__init__(**kwargs)

    def respond_reset_causality_indices(self, content):
        self.CA.reset_indices()
        reply = KQMLList('SUCCESS')
        return reply

    def respond_find_causal_path(self, content):
        """Response content to find-causal-path request"""

        source_arg = content.get('SOURCE')
        target_arg = content.get('TARGET')
        direction = content.gets('DIRECTION')

        if not source_arg or not target_arg:
            return self.make_failure('MISSING_MECHANISM')

        target_names = _get_kqml_names(target_arg)
        source_names = _get_kqml_names(source_arg)

        if not target_names or not source_names:
            return self.make_failure('NO_PATH_FOUND')

        source_name = source_names[0]
        target_name = target_names[0]

        target = {'id': target_name, 'pSite': ''}
        source = {'id': source_name, 'pSite': ''}

        result = self.CA.find_causality({'source': source, 'target': target, 'direction':direction})

        if not result:
            return self.make_failure('NO_PATH_FOUND')

        indra_json = [make_indra_json(result)]
        indra_stmts = stmts_from_json(indra_json)
        indra_cl_json = self.make_cljson(indra_stmts)

        reply = KQMLList('SUCCESS')
        reply.set('paths', indra_cl_json)

        # Send PC links to provenance tab
        self.send_provenance(result) # ['uri_str'])

        return reply

    def send_provenance(self, result):
        id1 = result['id1']
        mods1 = result['mods1']
        id2 = result['id2']
        mods2 = result['mods2']
        rel = result['rel']

        title = str(id1) +  ' ' + str(rel) + ' ' + str(id2)

        uri_str = result['uri_str']
        pc_url = 'http://www.pathwaycommons.org/pc2/get?' + uri_str + 'format=SBGN'
        html = '<a href= \'' + pc_url + '\' target= \'_blank\' > PC link</a>'
        msg = KQMLPerformative('tell')
        content = KQMLList('add-provenance')
        content.sets('html', html)
        pc_url_formatted = "http://www.pathwaycommons.org/pc2/get?" + uri_str + "format=SBGN"
        content.sets('pc', pc_url_formatted)
        content.sets('title', title)
        msg.set('content', content)
        self.send(msg)

    def respond_find_causality_target(self, content):
        """Response content to find-causality-target request"""
        target_arg = content.get('SOURCE')
        rel = content.gets('TYPE')

        if not target_arg:
            return self.make_failure('MISSING_MECHANISM')

        target_names = _get_kqml_names(target_arg)
        if not target_names:
            return self.make_failure('MISSING_MECHANISM')
        target_name = target_names[0]



        rel_map = {
            "phosphorylation": "phosphorylates",
            "dephosphorylation": "dephosphorylates",
            "activate": "upregulates-expression",
            "increase": "upregulates-expression",
            "inhibit": "downregulates-expression",
            "decrease": "downregulates-expression",
            "modulate": "modulates",
        }

        try:
            rel_verb = rel_map[rel]
        except:
            return self.make_failure('MISSING_MECHANISM')

        target = {'id': target_name, 'pSite': ' ', 'rel': rel_verb}
        result = self.CA.find_causality_targets(target)

        if not result:
            return self.make_failure('NO_PATH_FOUND')

        # Send PC links to provenance tab
        # Multiple interactions are sent separately
        for r in result:
            self.send_provenance(r) # r['uri_str'])

        indra_json = [make_indra_json(r) for r in result]
        indra_stmts = stmts_from_json(indra_json)
        indra_cl_json = self.make_cljson(indra_stmts)

        reply = KQMLList('SUCCESS')
        reply.set('paths', indra_cl_json)

        return reply

    def respond_find_causality_source(self, content):
        """Response content to find-qca-path request"""
        source_arg = content.get('TARGET')
        rel = content.gets('TYPE')

        if not source_arg:
            return self.make_failure('MISSING_MECHANISM')

        source_names = _get_kqml_names(source_arg)
        if not source_names:
            return self.make_failure('MISSING_MECHANISM')
        source_name = source_names[0]

        rel_map = {
            "phosphorylation": "is-phosphorylated-by",
            "dephosphorylation": "is-dephosphorylated-by",
            "activate": "expression-is-up",
            "increase": "upregulates-expression",
            "inhibit": "downregulates-expression",
            "decrease": "downregulates-expression",
            "modulate": "modulates",
        }

        try:
            rel_verb = rel_map[rel]
        except:
            return self.make_failure('MISSING_MECHANISM')

        source = {'id': source_name, 'pSite': ' ','rel': rel_verb}

        result = self.CA.find_causality_targets(source)

        if not result:
            return self.make_failure('NO_PATH_FOUND')

        # Multiple interactions are sent separately
        for r in result:
            self.send_provenance(r) # ['uri_str'])

        indra_json = [make_indra_json(r) for r in result]
        indra_stmts = stmts_from_json(indra_json)
        indra_cl_json = self.make_cljson(indra_stmts)

        reply = KQMLList('SUCCESS')
        reply.sets('paths', indra_cl_json)

        return reply

    def respond_dataset_correlated_entity(self, content):
        """Response content to find-database-correlated-entity request"""
        source_arg = content.get('SOURCE')
        if not source_arg:
            return self.make_failure('MISSING_MECHANISM')

        source_names = _get_kqml_names(source_arg)
        if not source_names:
            return self.make_failure('MISSING_MECHANISM')

        source_name = source_names[0]
        res = self.CA.find_next_correlation(source_name)
        if res == '':
            return self.make_failure('NO_PATH_FOUND')

        reply = KQMLList('SUCCESS')

        reply.sets('target', res['id2'])
        # TODO: should be float as original or converted to string?
        reply.sets('correlation', str(res['correlation']))
        reply.sets('explainable', res['explainable'])

        return reply


    def respond_find_common_upstreams(self, content):
        """Response content to find-common-upstreams request"""

        genes_arg = content.get('GENES')

        if not genes_arg:
            return self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(genes_arg)

        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        # gene_list = []
        # for gene_name in gene_names:
        #     gene_list.append(str(gene_name))

        result = self.CA.find_common_upstreams(gene_names)

        if not result:
            return self.make_failure('NO_UPSTREAM_FOUND')

        reply = KQMLList('SUCCESS')
        upstreams = _get_genes_cljson(result)
        reply.set('upstreams', upstreams)

        return reply

    def respond_find_mutation_significance(self, content):
        """Response content to find-mutation-significance request"""
        gene_arg = content.get('GENE')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        disease_arg = content.get('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_kqml_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = _sanitize_disase_name(disease_names[0])
        disease_abbr = self.CA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.CA.find_mutation_significance(gene_name, disease_abbr)

        if not result:
            return self.make_failure('MISSING_MECHANISM')

        reply = KQMLList('SUCCESS')
        reply.sets('mutsig', result)

        return reply

    def respond_find_mutex(self, content):
        """Response content to find-mutex request"""

        gene_arg = content.get('GENE')

        if not gene_arg:
            return self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        gene_name = gene_names[0]

        disease_arg = content.get('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')


        disease_names = _get_kqml_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = _sanitize_disase_name(disease_names[0])
        disease_abbr = self.CA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.CA.find_mutex(gene_name, disease_abbr)

        if not result:
            return self.make_failure('NO_MUTEX_GENES_FOUND')

        reply = KQMLList('SUCCESS')

        mutex = KQMLList()

        # Reorganize result
        for r in result:
            groups = KQMLList()
            groups.set('score', r['score'])

            genes = KQMLList()
            for gene in r['group']:
                genes.append(gene)

            groups.set('group', genes)
            mutex.append(groups)
            # mutex.append(groups)
        # TODO: cljson?
        reply.set('mutex', mutex)

        return reply



    def respond_find_cellular_location_from_names(self, content):
        """Response content to find-cellular-location-from-names request where genes are given as a list of names
        """

        gene_names = content.get('GENES')

        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        if isinstance(gene_names, str):
            return self.make_failure('INVALID_FORMAT')

        result = self.CA.find_most_likely_cellular_location(gene_names)

        if not result:
            return self.make_failure('NO_COMMON_CELLULAR_LOCATION_FOUND')

        reply = KQMLList('SUCCESS')

        components = _get_default_list_cljson(result)
        reply.set('components', components)

        return reply


    def respond_find_cellular_location(self, content):
        """Response content to find-cellular-location request"""
        genes_arg = content.get('GENES')

        if not genes_arg:
            return self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(genes_arg)

        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')

        # gene_list = []
        # for gene_name in gene_names:
        #     gene_list.append(str(gene_name))

        # result = self.CA.find_most_likely_cellular_location(gene_names)

        if not result:
            return self.make_failure('NO_COMMON_CELLULAR_LOCATION_FOUND')

        reply = KQMLList('SUCCESS')

        components = _get_default_list_cljson(result)
        reply.set('components', components)

        return reply

    def respond_find_mutation_frequency(self, content):
        """Response content to find-mutation-frequency request"""
        gene_arg = content.gets('GENE')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_term_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        disease_arg = content.gets('DISEASE')
        if not disease_arg:
            return self.make_failure('MISSING_MECHANISM')

        disease_names = _get_term_names(disease_arg)
        if not disease_names:
            return self.make_failure('INVALID_DISEASE')

        disease_name = _sanitize_disase_name(disease_names[0])
        disease_abbr = self.CA.get_tcga_abbr(disease_name)
        if disease_abbr is None:
            return self.make_failure('INVALID_DISEASE')

        result = self.CA.find_mutation_frequency(gene_name, disease_abbr)

        if not result:
            return self.make_failure('MISSING_MECHANISM')

        reply = KQMLList('SUCCESS')
        reply.sets('mutfreq', result)

        return reply

    def respond_find_gene_summary(self, content):
        gene_arg = content.get('GENE')

        if not gene_arg:
            self.make_failure('MISSING_MECHANISM')

        gene_names = _get_kqml_names(gene_arg)
        if not gene_names:
            return self.make_failure('MISSING_MECHANISM')
        gene_name = gene_names[0]

        result = self.CA.find_gene_summary(gene_name)

        reply = KQMLList('SUCCESS')
        reply.sets('geneSummary', result)

        return reply

def _get_kqml_names(kqmlList):
    """Given a kqml list returns the names of sublists in the list"""
    if not kqmlList:
        return None

    arr = kqmlList.data;
    if len(arr) == 0:
        return []

    if not isinstance(arr[0], KQMLList):
        arr = [kqmlList]

    res = list(map(lambda kl: kl.get('NAME').string_value(), arr))

    return res


def _sanitize_disase_name(name):
    """Given a disease name returns the sanitized version of it"""
    sanitized_name = name.replace("-", " ").lower()
    return sanitized_name


def _get_term_names(term_str):
    """Given an ekb-xml returns the names of genes in a list"""

    tp = TripsProcessor(term_str)
    terms = tp.tree.findall('TERM')
    if not terms:
        return None

    agent_names = []
    for term in terms:
        term_id = term.attrib['id']
        agent = tp._get_agent_by_id(term_id, None)

        if agent is not None:
            if isinstance(agent, list):
                for a in agent:
                    if a.name:
                        agent_names.append(a.name)
            else:
                agent_names.append(agent.name)

    if len(agent_names) == 0:
        return None

    return agent_names


def make_indra_json(causality):
    """Convert causality response to indra format
        Causality format is (id1, res1, pos1, id2,res2, pos2, rel)"""

    causality['rel'] = causality['rel'].upper()

    indra_relation_map = {
        "PHOSPHORYLATES": "Phosphorylation",
        "IS-PHOSPHORYLATED-BY": "Phosphorylation",
        "IS-DEPHOSPHORYLATED-BY": "Dephosphorylation",
        "UPREGULATES-EXPRESSION": "IncreaseAmount",
        "EXPRESSION-IS-UPREGULATED-BY": "IncreaseAmount",
        "DOWNREGULATES-EXPRESSION": "DecreaseAmount",
        "EXPRESSION-IS-DOWNREGULATED-BY": "DecreaseAmount"
    }

    rel_type = indra_relation_map[causality['rel']]

    s, t = ('2', '1') if 'IS' in causality['rel'] else ('1', '2')
    subj, obj = ('enz', 'sub') if 'PHOSPHO' in causality['rel'] else \
                ('subj', 'obj')

    # if "PHOSPHO" in causality['rel']:  # phosphorylation
    indra_json = {'type': rel_type,
                  subj: {'name': causality['id%s' % s],
                         'mods': causality['mods%s' % s]},
                  obj: {'name': causality['id%s' % t]},
                  'residue': causality['mods%s' % t][0]['residue'],
                  'position': causality['mods%s' % t][0]['position']}

    return indra_json

def _get_default_list_cljson(names):
    agents = list(map(lambda n: Agent(n), names))
    return Bioagent.make_cljson(agents)

def _get_genes_cljson(gene_names):
    agents = list(map(lambda n: _get_agent_from_gene_name(n), gene_names))
    return Bioagent.make_cljson(agents)

def _get_agent_from_gene_name(gene_name):
    db_refs = {}
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)


    if hgnc_id:
        db_refs['HGNC'] = hgnc_id

        uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
        if uniprot_id:
            db_refs['UP'] = uniprot_id

    agent = Agent(gene_name, db_refs=db_refs)
    return agent




if __name__ == "__main__":
    CausalityModule(argv=sys.argv[1:])
