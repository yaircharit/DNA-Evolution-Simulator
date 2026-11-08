from modules import config_loader
from modules.mutagenesis_simulator import mutagenesis_simulator

if __name__ == "__main__":
    simulator = mutagenesis_simulator()
    simulator.run_simulation()
    
