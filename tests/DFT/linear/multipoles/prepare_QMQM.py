import yaml
input_QMQM = yaml.load(open('waterQMQM_base.yaml', 'r').read())
output_QM = yaml.load(open('log-waterQM.yaml', 'r').read())
input_QMQM['dft']['external_potential'] = output_QM['Multipole coefficients']
f = open('waterQMQM.yaml', 'w')
f.write(yaml.dump(input_QMQM))
