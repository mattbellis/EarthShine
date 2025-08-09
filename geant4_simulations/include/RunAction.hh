// RunAction.hh
#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include <G4UserRunAction.hh>
#include <fstream>

class RunAction : public G4UserRunAction {
  public:
    RunAction();
    ~RunAction() override;
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    std::ofstream& GetStream() { return fOut; }

  private:
    std::ofstream fOut;
};

#endif

