import re

import pprzlink.message
import pprzlink.messages_xml_map

def parse_msg(line):
    return re.findall('^(\S*) (\S*) (\S*) (.*)$', line)[0] # ts, sender, msg_id, payload

class PprzBag:
    def __init__(self, filename):
        self.filename = filename

    def read_messages(self, msg_class, msg_ids):
        with open(self.filename) as fp:
            msgs = {msg_id:pprzlink.message.PprzMessage(msg_class, msg_id) for msg_id in msg_ids}
            for cnt, line in enumerate(fp):
                ts, sender, msg_id, payload = parse_msg(line)
                try:
                    msgs[msg_id].ivy_string_to_payload(payload)
                    yield float(ts), int(sender), msgs[msg_id]
                except KeyError: pass

# /home/poine/src/pprzlink/message_definitions/v1.0/messages.xml
class DataSet:
    def __init__(self, msg_class='telemetry', msg_ids=['GPS', 'ACTUATORS', 'IMU_GYRO'], msg_cbks=None):
        self.msg_class, self.msg_ids = msg_class, msg_ids
        def rec_gps(msg): return [float(msg.alt)/1000, float(msg.climb)/100]
        def rec_actuators(msg): return [float(_v) for _v in msg.values]
        def rec_imu_gyro(msg): return [float(_v) for _v in [msg.gp, msg.gq, msg.gr]]
        self.cbks = [rec_gps, rec_actuators, rec_imu_gyro] if msg_cbks is None else msg_cbks
        
    def readPprzBag(self, filename):
        self.tss = [[] for _m in self.msg_ids]
        self.vals = [[] for _m in self.msg_ids]
        msg_id_of_name = {msg_id:_i for _i, msg_id in enumerate(self.msg_ids)}
        for ts, sender, msg in PprzBag(filename).read_messages(self.msg_class, self.msg_ids):
            _id = msg_id_of_name[msg.name]
            self.tss[_id].append(ts)
            self.vals[_id].append(self.cbks[_id](msg))
    
    def save(self, filename):
        print('saving to {}'.format(filename))
        np.savez(filename, tss=np.asarray(self.tss), vals=np.asarray(self.vals))

    def load(self, filename):
        print('loading from {}'.format(filename))
        data =  np.load(filename, allow_pickle=True)
        self.tss, self.vals = data['tss'], data['vals']
